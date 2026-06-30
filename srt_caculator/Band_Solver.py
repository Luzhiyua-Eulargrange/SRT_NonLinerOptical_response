"""Band solver for the 1D two-component continuum Hamiltonian."""

from __future__ import annotations

from typing import Mapping

import numpy as np

try:
    from .config import normalize_params
except ImportError:  # Allows running this file directly from current directory.
    from config import normalize_params


def _basis_index(component: str, ell_index: int, num_g: int) -> int:
    if component == "A":
        return ell_index
    if component == "B":
        return num_g + ell_index
    raise ValueError("component must be 'A' or 'B'")


def build_hamiltonian(k_point: float, params: Mapping | None = None) -> np.ndarray:
    """
    Construct H(k) in the truncated plane-wave basis.
    Basis ordering is 'A, ell=-L..L' followed by 'B, ell=-L..L'.
    """
    p = normalize_params(params)
    ells = p["ell_list"]
    num_g = p["NG"]
    H = np.zeros((p["Nb"], p["Nb"]), dtype=np.complex128)

    b = p["b"]
    hbar = p["hbar"]
    kinetic_prefactor = hbar * hbar / (2.0 * p["m"])

    for i, ell in enumerate(ells):
        G = float(ell) * b
        row_A = _basis_index("A", i, num_g)
        row_B = _basis_index("B", i, num_g)

        for j, ell_p in enumerate(ells):
            col_A = _basis_index("A", j, num_g)
            col_B = _basis_index("B", j, num_g)

            if ell == ell_p:
                H[row_A, col_A] += kinetic_prefactor * (k_point + G - p["kappa"]) ** 2 + p["v0"]
                H[row_B, col_B] += -kinetic_prefactor * (k_point + G) ** 2 - p["v0"]

            if ell == ell_p + 1:
                H[row_A, col_A] += p["vA"]
                H[row_B, col_B] += -p["vB"]
                H[row_A, col_B] += p["w1"]
                H[row_B, col_A] += np.conjugate(p["w2"])

            if ell == ell_p - 1:
                H[row_A, col_A] += np.conjugate(p["vA"])
                H[row_B, col_B] += -np.conjugate(p["vB"])
                H[row_A, col_B] += p["w2"]
                H[row_B, col_A] += np.conjugate(p["w1"])

    return H


def diagonalize(k_point: float, params: Mapping | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted eigenvalues and eigenvectors of H(k)."""
    H = build_hamiltonian(float(k_point), params)
    energies, eigenvectors = np.linalg.eigh(H)
    return energies, eigenvectors


def solve_band(k_point: float, params: Mapping | None = None) -> tuple[np.ndarray, np.ndarray]:
    #Compatibility alias for diagonalizing one k point.
    #Forwarding diagonalize
    return diagonalize(k_point, params)


def smooth_eigenvector_phases(eigenvectors: np.ndarray, periodic: bool = True) -> np.ndarray:
    """
    Align adjacent eigenvector phases along a one-dimensional k grid.

    The operation changes only each band's U(1) phase. Eigenvectors stay
    normalized and orthogonal, and the corresponding band energies are
    unchanged. When ``periodic`` is true, the residual phase mismatch between
    the last and first k points is distributed uniformly over the whole grid
    instead of being left as a boundary jump.
    """
    vectors = np.array(eigenvectors, dtype=np.complex128, copy=True)
    if vectors.ndim != 3 or vectors.shape[1] != vectors.shape[2]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    for ik in range(1, vectors.shape[0]):
        for band in range(vectors.shape[2]):
            overlap = np.vdot(vectors[ik - 1, :, band], vectors[ik, :, band])
            if abs(overlap) > 0.0:
                vectors[ik, :, band] *= np.exp(-1j * np.angle(overlap))

    if periodic and vectors.shape[0] > 1:
        num_k = vectors.shape[0]
        phase_ramp_index = np.arange(num_k, dtype=float)
        for band in range(vectors.shape[2]):
            boundary_overlap = np.vdot(vectors[-1, :, band], vectors[0, :, band])
            if abs(boundary_overlap) > 0.0:
                boundary_phase = np.angle(boundary_overlap)
                phase_ramp = np.exp(1j * phase_ramp_index * boundary_phase / num_k)
                vectors[:, :, band] *= phase_ramp[:, None]
    return vectors


def phase_smoothing_diagnostics(eigenvectors: np.ndarray, include_boundary: bool = True) -> dict[str, float]:
    """
    Return overlap diagnostics for checking phase smoothing effectiveness.

    For periodic smoothing, the neighbor phases are expected to be small and
    spread across the whole grid instead of being concentrated at the FBZ
    boundary.
    """
    vectors = np.asarray(eigenvectors, dtype=np.complex128)
    if vectors.ndim != 3 or vectors.shape[1] != vectors.shape[2]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")
    if vectors.shape[0] < 2:
        raise ValueError("at least two k points are required")

    overlaps = []
    phases = []
    for ik in range(1, vectors.shape[0]):
        for band in range(vectors.shape[2]):
            overlap = np.vdot(vectors[ik - 1, :, band], vectors[ik, :, band])
            overlaps.append(abs(overlap))
            phases.append(abs(np.angle(overlap)) if abs(overlap) > 0.0 else np.nan)

    overlap_array = np.asarray(overlaps, dtype=float)
    phase_array = np.asarray(phases, dtype=float)
    finite_phases = phase_array[np.isfinite(phase_array)]
    diagnostics = {
        "min_adjacent_overlap_abs": float(np.min(overlap_array)),
        "mean_adjacent_overlap_abs": float(np.mean(overlap_array)),
        "max_adjacent_phase_abs": float(np.max(finite_phases)) if finite_phases.size else float("nan"),
        "mean_adjacent_phase_abs": float(np.mean(finite_phases)) if finite_phases.size else float("nan"),
    }
    if include_boundary:
        boundary_overlaps = []
        boundary_phases = []
        for band in range(vectors.shape[2]):
            overlap = np.vdot(vectors[-1, :, band], vectors[0, :, band])
            boundary_overlaps.append(abs(overlap))
            boundary_phases.append(abs(np.angle(overlap)) if abs(overlap) > 0.0 else np.nan)
        boundary_overlap_array = np.asarray(boundary_overlaps, dtype=float)
        boundary_phase_array = np.asarray(boundary_phases, dtype=float)
        finite_boundary_phases = boundary_phase_array[np.isfinite(boundary_phase_array)]
        diagnostics.update({
            "min_boundary_overlap_abs": float(np.min(boundary_overlap_array)),
            "mean_boundary_overlap_abs": float(np.mean(boundary_overlap_array)),
            "max_boundary_phase_abs": (
                float(np.max(finite_boundary_phases)) if finite_boundary_phases.size else float("nan")
            ),
            "mean_boundary_phase_abs": (
                float(np.mean(finite_boundary_phases)) if finite_boundary_phases.size else float("nan")
            ),
        })
    return diagnostics


def solve_bands(
    k_grid: np.ndarray,
    params: Mapping | None = None,
    smooth_phases: bool = True,
    periodic_smooth: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Diagonalize H(k) on a k grid.
    Returns energies with shape (Nk, Nb) and eigenvectors with shape
    (Nk, Nb, Nb).  For each k, eigenvectors are stored column-wise.
    By default, eigenvector phases are smoothed along the k grid before return,
    including periodic phase-ramp smoothing across the FBZ boundary.
    """
    p = normalize_params(params)
    k_grid = np.asarray(k_grid, dtype=float)
    energies = np.empty((k_grid.size, p["Nb"]), dtype=float)
    eigenvectors = np.empty((k_grid.size, p["Nb"], p["Nb"]), dtype=np.complex128)

    for ik, k_point in enumerate(k_grid):
        energies[ik], eigenvectors[ik] = diagonalize(float(k_point), p)

    if smooth_phases:
        eigenvectors = smooth_eigenvector_phases(eigenvectors, periodic=periodic_smooth)

    return energies, eigenvectors


def save_band_results(
    filename: str,
    k_grid: np.ndarray,
    k_weight: float,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
) -> None:
    """Save band calculation inputs and outputs to an npz file."""
    np.savez(
        filename,
        k_grid=np.asarray(k_grid, dtype=float),
        k_weight=float(k_weight),
        energies=np.asarray(energies, dtype=float),
        eigenvectors=np.asarray(eigenvectors, dtype=np.complex128),
    )


def load_band_results(filename: str) -> dict[str, np.ndarray | float]:
    """Load band calculation data saved by save_band_results."""
    data = np.load(filename)
    return {
        "k_grid": data["k_grid"],
        "k_weight": float(data["k_weight"]),
        "energies": data["energies"],
        "eigenvectors": data["eigenvectors"],
    }


def is_hermitian(matrix: np.ndarray, atol: float = 1e-10) -> bool:
    """Return True if a matrix is Hermitian within numerical tolerance."""
    return bool(np.allclose(matrix, matrix.conj().T, atol=atol))
