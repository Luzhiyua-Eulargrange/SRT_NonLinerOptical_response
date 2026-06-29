"""
Independent diagnostics for comparing band-velocity formulas.

This module is intentionally not imported by the main RDM workflow. The
production velocity-gauge solver uses the analytic band-basis transform
U^dagger (1/hbar dH/dk) U. The finite-difference Eq. (28) construction below is
kept only to study gauge-smoothing and Berry-connection numerical issues.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

try:
    from .Band_Solver import solve_bands
    from .Geometry import berry_connection, velocity_matrix
    from .config import make_k_grid, normalize_params
except ImportError:
    from Band_Solver import solve_bands
    from Geometry import berry_connection, velocity_matrix
    from config import make_k_grid, normalize_params


def velocity_from_eq28_at_indices(
    k_grid: np.ndarray,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
    params: Mapping | None = None,
    indices: tuple[int, ...] = (0,),
) -> tuple[np.ndarray, np.ndarray]:
    """Compute Eq. (28) velocity matrices only at selected k indices."""
    p = normalize_params(params)
    k_grid = np.asarray(k_grid, dtype=float)
    energies = np.asarray(energies, dtype=float)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)

    if k_grid.ndim != 1 or k_grid.size < 3:
        raise ValueError("k_grid must be one-dimensional with at least three points")
    if energies.ndim != 2 or energies.shape[0] != k_grid.size:
        raise ValueError("energies must have shape (Nk, Nb)")
    if vectors.shape != (k_grid.size, energies.shape[1], energies.shape[1]):
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    steps = np.diff(k_grid)
    if not np.allclose(steps, steps[0]):
        raise ValueError("Eq. (28) diagnostic requires a uniform k grid")

    nk, nb = energies.shape
    dk = float(steps[0])
    xi_grid = berry_connection(k_grid, vectors, smooth_gauge=False)
    diag = np.arange(nb)
    normalized_indices: list[int] = []
    velocities: list[np.ndarray] = []

    for index in indices:
        ik = int(index)
        if ik < 0:
            ik += nk
        if ik < 0 or ik >= nk:
            raise IndexError(f"k index {index} is out of range for Nk={nk}")

        im = (ik - 1) % nk
        ip = (ik + 1) % nk
        d_energy = (energies[ip] - energies[im]) / (2.0 * dk)
        xi = xi_grid[ik]

        energy_diff = energies[ik, None, :] - energies[ik, :, None]
        velocity = -1j * energy_diff * xi
        velocity[diag, diag] = d_energy
        normalized_indices.append(ik)
        velocities.append(velocity / p["hbar"])

    return np.asarray(normalized_indices, dtype=int), np.asarray(velocities, dtype=np.complex128)


def velocity_from_basis_transform_at_indices(
    k_grid: np.ndarray,
    eigenvectors: np.ndarray,
    params: Mapping | None = None,
    indices: tuple[int, ...] = (0,),
) -> tuple[np.ndarray, np.ndarray]:
    """Compute U^dagger (1/hbar dH/dk) U only at selected k indices."""
    p = normalize_params(params)
    k_grid = np.asarray(k_grid, dtype=float)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)
    if k_grid.ndim != 1:
        raise ValueError("k_grid must be one-dimensional")
    if vectors.ndim != 3 or vectors.shape[0] != k_grid.size or vectors.shape[1] != vectors.shape[2]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    normalized_indices: list[int] = []
    velocities: list[np.ndarray] = []
    for index in indices:
        ik = int(index)
        if ik < 0:
            ik += k_grid.size
        if ik < 0 or ik >= k_grid.size:
            raise IndexError(f"k index {index} is out of range for Nk={k_grid.size}")

        v_plane = velocity_matrix(float(k_grid[ik]), p, basis="plane_wave")
        u = vectors[ik]
        normalized_indices.append(ik)
        velocities.append(u.conj().T @ v_plane @ u)

    return np.asarray(normalized_indices, dtype=int), np.asarray(velocities, dtype=np.complex128)


def compare_velocity_formulas(
    k_grid: np.ndarray,
    params: Mapping | None,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
    indices: tuple[int, ...] = (0,),
) -> dict[str, float]:
    """Compare Eq. (28) and basis-transform velocities at selected k indices."""
    selected_a, velocity_eq28 = velocity_from_eq28_at_indices(
        k_grid,
        energies,
        eigenvectors,
        params=params,
        indices=indices,
    )
    selected_b, velocity_transform = velocity_from_basis_transform_at_indices(
        k_grid,
        eigenvectors,
        params=params,
        indices=indices,
    )
    if not np.array_equal(selected_a, selected_b):
        raise RuntimeError("diagnostic index normalization mismatch")

    difference = velocity_eq28 - velocity_transform
    max_reference = float(np.max(np.abs(velocity_transform)))
    max_difference = float(np.max(np.abs(difference)))
    return {
        "num_checked": float(selected_a.size),
        "max_abs_difference": max_difference,
        "rms_difference": float(np.sqrt(np.mean(np.abs(difference) ** 2))),
        "max_reference_abs": max_reference,
        "relative_max_difference": max_difference / max(max_reference, 1e-30),
    }


def run_default_diagnostic() -> dict[str, float]:
    """Run a standalone velocity-matrix diagnostic with conservative defaults."""
    params = normalize_params({
        "L": 2,
        "num_k_rdm": 101,
    })
    k_grid, _ = make_k_grid(params)
    energies, eigenvectors = solve_bands(k_grid, params)
    indices = (1, len(k_grid) // 4, len(k_grid) // 2, 3 * len(k_grid) // 4, -2)
    result = compare_velocity_formulas(
        k_grid,
        params,
        energies,
        eigenvectors,
        indices=indices,
    )
    result["num_k"] = float(k_grid.size)
    return result


def main() -> None:
    """Entry point for independent diagnostics."""
    result = run_default_diagnostic()
    print("Velocity matrix diagnostic")
    print("  method A: U^dagger (1/hbar dH/dk) U")
    print("  method B: Eq. (28) finite-difference Berry connection")
    print("  num_k:", int(result["num_k"]))
    print("  checked k points:", int(result["num_checked"]))
    print("  max |method B - method A|:", result["max_abs_difference"])
    print("  rms |method B - method A|:", result["rms_difference"])
    print("  max |method A|:", result["max_reference_abs"])
    print("  relative max difference:", result["relative_max_difference"])
    print("  note: this diagnostic is not used by main.py or RDM evolution.")


if __name__ == "__main__":
    main()
