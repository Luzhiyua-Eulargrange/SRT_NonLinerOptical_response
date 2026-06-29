"""Unit-cell real-space density reconstructed from band-basis RDM data."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np

try:
    from .config import normalize_params
except ImportError:
    from config import normalize_params


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


def make_real_space_grid(params: Mapping | None, num_x: int = 201) -> tuple[np.ndarray, float]:
    """Return an endpoint-excluding grid over one unit cell and its dx."""
    p = normalize_params(params)
    if int(num_x) < 2:
        raise ValueError("num_x must be at least 2")
    x_grid = np.linspace(0.0, p["a0"], int(num_x), endpoint=False)
    return x_grid, p["a0"] / int(num_x)


def band_to_plane_wave_density(
    rho_band_trajectory: np.ndarray,
    eigenvectors: np.ndarray,
) -> np.ndarray:
    """
    Transform rho from band basis to the plane-wave basis.

    ``eigenvectors[k]`` stores plane-wave basis vectors column-wise, so
    ``rho_pw(k,t) = U(k) rho_band(k,t) U(k)^dagger``.
    """
    rho_band = np.asarray(rho_band_trajectory, dtype=np.complex128)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)
    if rho_band.ndim != 4 or rho_band.shape[-1] != rho_band.shape[-2]:
        raise ValueError("rho_band_trajectory must have shape (Nt, Nk, Nb, Nb)")
    if vectors.shape != rho_band.shape[1:]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    return np.einsum("kis,tksu,kju->tkij", vectors, rho_band, vectors.conj(), optimize=True)


def unit_cell_density_from_band_rdm(
    params: Mapping | None,
    k_grid: np.ndarray,
    k_weight: float,
    rho_trajectory: np.ndarray,
    eigenvectors: np.ndarray,
    x_grid: np.ndarray | None = None,
    num_x: int = 201,
) -> dict:
    """
    Reconstruct component-resolved density inside one unit cell.

    The returned density is periodic over one cell. It assumes the RDM is block
    diagonal in crystal momentum, as is appropriate for a spatially uniform
    external field and a translation-invariant initial state.
    """
    k_grid = np.asarray(k_grid, dtype=float)
    rho = np.asarray(rho_trajectory, dtype=np.complex128)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)
    if rho.ndim != 4 or rho.shape[-1] != rho.shape[-2]:
        raise ValueError("rho_trajectory must have shape (Nt, Nk, Nb, Nb)")
    if k_grid.ndim != 1 or k_grid.size != rho.shape[1]:
        raise ValueError("k_grid size must match rho_trajectory")
    if vectors.shape != rho.shape[1:]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    p = _params_for_basis(params, rho.shape[-1])
    if x_grid is None:
        x_grid, dx = make_real_space_grid(p, num_x=num_x)
    else:
        x_grid = np.asarray(x_grid, dtype=float)
        if x_grid.ndim != 1 or x_grid.size < 2:
            raise ValueError("x_grid must be one-dimensional with at least two points")
        dx = p["a0"] / x_grid.size

    reciprocal_vectors = p["ell_list"].astype(float) * p["b"]
    phase = np.exp(1j * (reciprocal_vectors[:, None] - reciprocal_vectors[None, :])[None, :, :] * x_grid[:, None, None])
    rho_pw = band_to_plane_wave_density(rho, vectors)

    nt = rho.shape[0]
    num_g = p["NG"]
    density_components = np.empty((nt, 2, x_grid.size), dtype=float)
    for component_index, start in enumerate((0, num_g)):
        block = rho_pw[:, :, start:start + num_g, start:start + num_g]
        density_components[:, component_index] = (
            float(k_weight) * np.einsum("tkij,xij->tx", block, phase, optimize=True)
        ).real

    density_total = np.sum(density_components, axis=1)
    charge_density = -p["e_charge"] * density_total
    component_charge_density = -p["e_charge"] * density_components
    particles_per_cell = dx * np.sum(density_total, axis=1)
    charge_per_cell = dx * np.sum(charge_density, axis=1)

    return {
        "x_grid": x_grid,
        "dx": float(dx),
        "k_grid": k_grid,
        "k_weight": float(k_weight),
        "time_grid": None,
        "density_components": density_components,
        "density_total": density_total,
        "charge_density_components": component_charge_density,
        "charge_density_total": charge_density,
        "particles_per_cell": particles_per_cell,
        "charge_per_cell": charge_per_cell,
    }


def density_result_from_rdm_result(
    params: Mapping | None,
    rdm_result: Mapping,
    x_grid: np.ndarray | None = None,
    num_x: int = 201,
) -> dict:
    """Return a unit-cell density result dictionary from an RDM result mapping."""
    required_keys = ("k_grid", "k_weight", "time_grid", "rho_trajectory", "eigenvectors")
    missing = [key for key in required_keys if key not in rdm_result]
    if missing:
        raise KeyError(f"rdm_result is missing keys: {missing}")

    result = unit_cell_density_from_band_rdm(
        params,
        rdm_result["k_grid"],
        float(rdm_result["k_weight"]),
        rdm_result["rho_trajectory"],
        rdm_result["eigenvectors"],
        x_grid=x_grid,
        num_x=num_x,
    )
    result["time_grid"] = np.asarray(rdm_result["time_grid"], dtype=float)
    return result


def density_convergence(reference: Mapping, candidate: Mapping) -> dict[str, float]:
    """Compare two density results on the same x and time grids."""
    for key in ("x_grid", "time_grid"):
        if not np.allclose(reference[key], candidate[key], rtol=0.0, atol=1e-12):
            raise ValueError(f"density results use different {key}")

    ref_total = np.asarray(reference["density_total"], dtype=float)
    cand_total = np.asarray(candidate["density_total"], dtype=float)
    ref_components = np.asarray(reference["density_components"], dtype=float)
    cand_components = np.asarray(candidate["density_components"], dtype=float)
    ref_charge = np.asarray(reference["charge_per_cell"], dtype=float)
    cand_charge = np.asarray(candidate["charge_per_cell"], dtype=float)

    total_scale = max(float(np.max(np.abs(ref_total))), 1e-12)
    component_scale = max(float(np.max(np.abs(ref_components))), 1e-12)
    charge_scale = max(float(np.max(np.abs(ref_charge))), 1e-12)
    return {
        "density_total_relative_max_difference": float(np.max(np.abs(cand_total - ref_total)) / total_scale),
        "density_total_max_abs_difference": float(np.max(np.abs(cand_total - ref_total))),
        "density_component_relative_max_difference": float(np.max(np.abs(cand_components - ref_components)) / component_scale),
        "density_component_max_abs_difference": float(np.max(np.abs(cand_components - ref_components))),
        "charge_per_cell_relative_max_difference": float(np.max(np.abs(cand_charge - ref_charge)) / charge_scale),
        "charge_per_cell_max_abs_difference": float(np.max(np.abs(cand_charge - ref_charge))),
    }


def save_density_result(filename: str, density_result: Mapping) -> None:
    """Save unit-cell density data to an npz file."""
    np.savez(
        filename,
        x_grid=np.asarray(density_result["x_grid"], dtype=float),
        dx=float(density_result["dx"]),
        k_grid=np.asarray(density_result["k_grid"], dtype=float),
        k_weight=float(density_result["k_weight"]),
        time_grid=np.asarray(density_result["time_grid"], dtype=float),
        density_components=np.asarray(density_result["density_components"], dtype=float),
        density_total=np.asarray(density_result["density_total"], dtype=float),
        charge_density_components=np.asarray(density_result["charge_density_components"], dtype=float),
        charge_density_total=np.asarray(density_result["charge_density_total"], dtype=float),
        particles_per_cell=np.asarray(density_result["particles_per_cell"], dtype=float),
        charge_per_cell=np.asarray(density_result["charge_per_cell"], dtype=float),
    )


def load_rdm_result(filename: str | Path) -> dict:
    """Load an RDM npz file into a plain dictionary."""
    data = np.load(filename)
    return {key: data[key] for key in data.files}


def main() -> None:
    """Build density data from the default length-gauge RDM result file."""
    base_dir = Path(__file__).resolve().parent
    input_path = base_dir / "length_rdm_result.npz"
    output_path = base_dir / "length_density_result.npz"
    if not input_path.exists():
        raise FileNotFoundError(
            f"Could not find {input_path}. Run python main.py first to generate the length-gauge RDM result."
        )

    rdm_result = load_rdm_result(input_path)
    density_result = density_result_from_rdm_result(None, rdm_result)
    save_density_result(output_path, density_result)
    print(f"Saved density result to {output_path}")
    print("  density_total shape:", density_result["density_total"].shape)
    print("  density_components shape:", density_result["density_components"].shape)
    print("  particles_per_cell drift:", float(np.max(np.abs(
        density_result["particles_per_cell"] - density_result["particles_per_cell"][0]
    ))))


if __name__ == "__main__":
    main()
