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


def solve_bands(k_grid: np.ndarray, params: Mapping | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Diagonalize H(k) on a k grid.
    Returns energies with shape (Nk, Nb) and eigenvectors with shape
    (Nk, Nb, Nb).  For each k, eigenvectors are stored column-wise.
    """
    p = normalize_params(params)
    k_grid = np.asarray(k_grid, dtype=float)
    energies = np.empty((k_grid.size, p["Nb"]), dtype=float)
    eigenvectors = np.empty((k_grid.size, p["Nb"], p["Nb"]), dtype=np.complex128)

    for ik, k_point in enumerate(k_grid):
        energies[ik], eigenvectors[ik] = diagonalize(float(k_point), p)

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
