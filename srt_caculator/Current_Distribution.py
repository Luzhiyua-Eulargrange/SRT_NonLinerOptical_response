"""Unit-cell real-space current reconstructed from band-basis RDM data.

The implementation follows the local current operator derived by minimal
coupling in the plane-wave basis. It does not use perturbative energy
denominators or k-derivatives of eigenvectors.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping

import numpy as np

try:
    from .Density import band_to_plane_wave_density, load_rdm_result, make_real_space_grid
    from .config import normalize_params
except ImportError:
    from Density import band_to_plane_wave_density, load_rdm_result, make_real_space_grid
    from config import normalize_params


VectorPotential = None | float | np.ndarray | Callable[[float], float]


def _params_for_basis(params: Mapping | None, num_bands: int) -> dict:
    """Return validated parameters whose basis size matches the RDM data."""
    if num_bands % 2 != 0:
        raise ValueError("num_bands must be even for the two-component basis")

    num_g = num_bands // 2
    if num_g % 2 != 1:
        raise ValueError("number of reciprocal vectors must be odd")

    inferred_l = (num_g - 1) // 2
    if params is None:
        return normalize_params({"L": inferred_l})

    p = normalize_params(params)
    if p["Nb"] != num_bands:
        raise ValueError(f"params imply Nb={p['Nb']}, but RDM data has Nb={num_bands}")
    return p


def vector_potential_on_grid(
    time_grid: np.ndarray | None,
    vector_potential: VectorPotential = None,
) -> np.ndarray:
    """Return A(t) sampled on time_grid.

    ``None`` means the length-gauge/default convention A(t)=0. A scalar gives a
    constant value, an array must match time_grid, and a callable is sampled at
    every time point.
    """
    if time_grid is None:
        if vector_potential is None:
            return np.zeros(1, dtype=float)
        if np.isscalar(vector_potential):
            return np.asarray([float(vector_potential)], dtype=float)
        values = np.asarray(vector_potential, dtype=float)
        if values.ndim != 1:
            raise ValueError("vector_potential array must be one-dimensional")
        return values

    times = np.asarray(time_grid, dtype=float)
    if times.ndim != 1 or times.size < 1:
        raise ValueError("time_grid must be a one-dimensional array")

    if vector_potential is None:
        return np.zeros(times.size, dtype=float)
    if np.isscalar(vector_potential):
        return np.full(times.size, float(vector_potential), dtype=float)
    if callable(vector_potential):
        return np.asarray([float(vector_potential(float(t))) for t in times], dtype=float)

    values = np.asarray(vector_potential, dtype=float)
    if values.shape != times.shape:
        raise ValueError("vector_potential array must have the same shape as time_grid")
    return values


def _validate_current_inputs(
    k_grid: np.ndarray,
    rho_trajectory: np.ndarray,
    eigenvectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k_grid = np.asarray(k_grid, dtype=float)
    rho = np.asarray(rho_trajectory, dtype=np.complex128)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)

    if k_grid.ndim != 1:
        raise ValueError("k_grid must be one-dimensional")
    if rho.ndim != 4 or rho.shape[-1] != rho.shape[-2]:
        raise ValueError("rho_trajectory must have shape (Nt, Nk, Nb, Nb)")
    if rho.shape[1] != k_grid.size:
        raise ValueError("k_grid size must match rho_trajectory")
    if vectors.shape != rho.shape[1:]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")
    return k_grid, rho, vectors


def unit_cell_current_from_band_rdm(
    params: Mapping | None,
    k_grid: np.ndarray,
    k_weight: float,
    rho_trajectory: np.ndarray,
    eigenvectors: np.ndarray,
    time_grid: np.ndarray | None = None,
    x_grid: np.ndarray | None = None,
    num_x: int = 201,
    vector_potential: VectorPotential = None,
) -> dict:
    """
    Reconstruct component-resolved charge-current density inside one unit cell.

    The plane-wave convention is the same as ``Density.py``:

        u_{ks,alpha}(x) = sum_G U_{alpha G,s}(k) exp(i G x)

    and the Brillouin-zone integral is discretized with ``k_weight = dk/(2*pi)``.
    The returned current is the charge current for electron charge ``-e`` with
    ``e = params["e_charge"] > 0``.
    """
    k_grid, rho, vectors = _validate_current_inputs(k_grid, rho_trajectory, eigenvectors)
    p = _params_for_basis(params, rho.shape[-1])

    if time_grid is not None:
        times = np.asarray(time_grid, dtype=float)
        if times.ndim != 1 or times.size != rho.shape[0]:
            raise ValueError("time_grid must have shape (Nt,)")
    else:
        times = None

    if x_grid is None:
        x_grid, dx = make_real_space_grid(p, num_x=num_x)
    else:
        x_grid = np.asarray(x_grid, dtype=float)
        if x_grid.ndim != 1 or x_grid.size < 2:
            raise ValueError("x_grid must be one-dimensional with at least two points")
        dx = p["a0"] / x_grid.size

    a_values = vector_potential_on_grid(times, vector_potential)
    if a_values.size != rho.shape[0]:
        if a_values.size == 1:
            a_values = np.full(rho.shape[0], float(a_values[0]), dtype=float)
        else:
            raise ValueError("vector_potential does not match the number of time steps")

    rho_pw = band_to_plane_wave_density(rho, vectors)

    num_g = p["NG"]
    reciprocal_vectors = p["ell_list"].astype(float) * p["b"]
    phase = np.exp(
        1j
        * (reciprocal_vectors[:, None] - reciprocal_vectors[None, :])[None, :, :]
        * x_grid[:, None, None]
    )

    momentum_sum = (
        2.0 * k_grid[:, None, None]
        + reciprocal_vectors[None, :, None]
        + reciprocal_vectors[None, None, :]
    )
    prefactor = p["e_charge"] * p["hbar"] / (2.0 * p["m"])
    coeff_a = prefactor * (-momentum_sum + 2.0 * p["kappa"])
    coeff_b = prefactor * momentum_sum

    current_para = np.empty((rho.shape[0], 2, x_grid.size), dtype=float)
    density_components = np.empty_like(current_para)
    for component_index, start in enumerate((0, num_g)):
        block = rho_pw[:, :, start:start + num_g, start:start + num_g]
        coeff = coeff_a if component_index == 0 else coeff_b
        current_para[:, component_index] = (
            float(k_weight) * np.einsum("tkij,kij,xij->tx", block, coeff, phase, optimize=True)
        ).real
        density_components[:, component_index] = (
            float(k_weight) * np.einsum("tkij,xij->tx", block, phase, optimize=True)
        ).real

    diamagnetic_coeff = np.empty((rho.shape[0], 2), dtype=float)
    diamagnetic_coeff[:, 0] = -(p["e_charge"] ** 2) * a_values / p["m"]
    diamagnetic_coeff[:, 1] = +(p["e_charge"] ** 2) * a_values / p["m"]
    current_dia = diamagnetic_coeff[:, :, None] * density_components

    current_components = current_para + current_dia
    current_total = np.sum(current_components, axis=1)
    paramagnetic_current_total = np.sum(current_para, axis=1)
    diamagnetic_current_total = np.sum(current_dia, axis=1)
    current_per_cell = float(dx) * np.sum(current_total, axis=1)
    average_current = current_per_cell / p["a0"]

    return {
        "x_grid": x_grid,
        "dx": float(dx),
        "k_grid": k_grid,
        "k_weight": float(k_weight),
        "time_grid": None if times is None else times,
        "vector_potential": a_values,
        "density_components": density_components,
        "paramagnetic_current_components": current_para,
        "diamagnetic_current_components": current_dia,
        "current_components": current_components,
        "current_total": current_total,
        "paramagnetic_current_total": paramagnetic_current_total,
        "diamagnetic_current_total": diamagnetic_current_total,
        "current_per_cell": current_per_cell,
        "average_current": average_current,
    }


def current_result_from_rdm_result(
    params: Mapping | None,
    rdm_result: Mapping,
    x_grid: np.ndarray | None = None,
    num_x: int = 201,
    vector_potential: VectorPotential = None,
) -> dict:
    """Return a unit-cell current result dictionary from an RDM result mapping."""
    required_keys = ("k_grid", "k_weight", "rho_trajectory", "eigenvectors")
    missing = [key for key in required_keys if key not in rdm_result]
    if missing:
        raise KeyError(f"rdm_result is missing keys: {missing}")

    time_grid = rdm_result["time_grid"] if "time_grid" in rdm_result else None
    return unit_cell_current_from_band_rdm(
        params,
        rdm_result["k_grid"],
        float(rdm_result["k_weight"]),
        rdm_result["rho_trajectory"],
        rdm_result["eigenvectors"],
        time_grid=time_grid,
        x_grid=x_grid,
        num_x=num_x,
        vector_potential=vector_potential,
    )


def _relative_max_difference(reference: np.ndarray, candidate: np.ndarray, floor: float = 1e-12) -> float:
    """Return max-norm relative difference between two arrays."""
    reference = np.asarray(reference)
    candidate = np.asarray(candidate)
    if reference.shape != candidate.shape:
        raise ValueError("arrays must have the same shape")
    scale = max(float(np.max(np.abs(reference))), float(floor))
    return float(np.max(np.abs(candidate - reference)) / scale)


def current_convergence(reference: Mapping, candidate: Mapping) -> dict[str, float]:
    """Compare two current-distribution results on the same x and time grids."""
    for key in ("x_grid", "time_grid"):
        if reference[key] is None or candidate[key] is None:
            if reference[key] is not candidate[key]:
                raise ValueError(f"current results use different {key}")
        elif not np.allclose(reference[key], candidate[key], rtol=0.0, atol=1e-12):
            raise ValueError(f"current results use different {key}")

    for key in ("vector_potential",):
        if not np.allclose(reference[key], candidate[key], rtol=0.0, atol=1e-12):
            raise ValueError(f"current results use different {key}")

    ref_total = np.asarray(reference["current_total"], dtype=float)
    cand_total = np.asarray(candidate["current_total"], dtype=float)
    ref_components = np.asarray(reference["current_components"], dtype=float)
    cand_components = np.asarray(candidate["current_components"], dtype=float)
    ref_para = np.asarray(reference["paramagnetic_current_components"], dtype=float)
    cand_para = np.asarray(candidate["paramagnetic_current_components"], dtype=float)
    ref_dia = np.asarray(reference["diamagnetic_current_components"], dtype=float)
    cand_dia = np.asarray(candidate["diamagnetic_current_components"], dtype=float)
    ref_per_cell = np.asarray(reference["current_per_cell"], dtype=float)
    cand_per_cell = np.asarray(candidate["current_per_cell"], dtype=float)
    ref_average = np.asarray(reference["average_current"], dtype=float)
    cand_average = np.asarray(candidate["average_current"], dtype=float)

    return {
        "current_total_relative_max_difference": _relative_max_difference(ref_total, cand_total),
        "current_total_max_abs_difference": float(np.max(np.abs(cand_total - ref_total))),
        "current_component_relative_max_difference": _relative_max_difference(ref_components, cand_components),
        "current_component_max_abs_difference": float(np.max(np.abs(cand_components - ref_components))),
        "paramagnetic_component_relative_max_difference": _relative_max_difference(ref_para, cand_para),
        "paramagnetic_component_max_abs_difference": float(np.max(np.abs(cand_para - ref_para))),
        "diamagnetic_component_relative_max_difference": _relative_max_difference(ref_dia, cand_dia),
        "diamagnetic_component_max_abs_difference": float(np.max(np.abs(cand_dia - ref_dia))),
        "current_per_cell_relative_max_difference": _relative_max_difference(ref_per_cell, cand_per_cell),
        "current_per_cell_max_abs_difference": float(np.max(np.abs(cand_per_cell - ref_per_cell))),
        "average_current_relative_max_difference": _relative_max_difference(ref_average, cand_average),
        "average_current_max_abs_difference": float(np.max(np.abs(cand_average - ref_average))),
    }


def save_current_result(filename: str | Path, current_result: Mapping) -> None:
    """Save unit-cell current data to an npz file."""
    time_grid = current_result["time_grid"]
    np.savez(
        filename,
        x_grid=np.asarray(current_result["x_grid"], dtype=float),
        dx=float(current_result["dx"]),
        k_grid=np.asarray(current_result["k_grid"], dtype=float),
        k_weight=float(current_result["k_weight"]),
        time_grid=np.asarray([] if time_grid is None else time_grid, dtype=float),
        vector_potential=np.asarray(current_result["vector_potential"], dtype=float),
        density_components=np.asarray(current_result["density_components"], dtype=float),
        paramagnetic_current_components=np.asarray(
            current_result["paramagnetic_current_components"], dtype=float
        ),
        diamagnetic_current_components=np.asarray(
            current_result["diamagnetic_current_components"], dtype=float
        ),
        current_components=np.asarray(current_result["current_components"], dtype=float),
        current_total=np.asarray(current_result["current_total"], dtype=float),
        paramagnetic_current_total=np.asarray(
            current_result["paramagnetic_current_total"], dtype=float
        ),
        diamagnetic_current_total=np.asarray(
            current_result["diamagnetic_current_total"], dtype=float
        ),
        current_per_cell=np.asarray(current_result["current_per_cell"], dtype=float),
        average_current=np.asarray(current_result["average_current"], dtype=float),
    )


def main() -> None:
    """Build current data from the default length-gauge RDM result file."""
    base_dir = Path(__file__).resolve().parent
    input_path = base_dir / "length_rdm_result.npz"
    output_path = base_dir / "length_current_distribution.npz"
    if not input_path.exists():
        raise FileNotFoundError(
            f"Could not find {input_path}. Run python main.py first to generate the length-gauge RDM result."
        )

    rdm_result = load_rdm_result(input_path)
    current_result = current_result_from_rdm_result(None, rdm_result)
    save_current_result(output_path, current_result)
    print(f"Saved current distribution to {output_path}")
    print("  current_total shape:", current_result["current_total"].shape)
    print("  current_components shape:", current_result["current_components"].shape)
    print("  max |average_current|:", float(np.max(np.abs(current_result["average_current"]))))


if __name__ == "__main__":
    main()
