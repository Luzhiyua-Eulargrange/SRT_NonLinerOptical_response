"""Berry-connection geometry helpers for the 1D continuum model."""

from __future__ import annotations

import numpy as np

try:
    from .Band_Solver import smooth_eigenvector_phases
except ImportError:
    from Band_Solver import smooth_eigenvector_phases


def berry_connection(
    k_grid: np.ndarray,
    eigenvectors: np.ndarray,
    smooth_gauge: bool = True,
) -> np.ndarray:
    # Compute xi_ss'(k) = i <u_ks | d_k u_ks'> by finite differences.
    k_grid = np.asarray(k_grid, dtype=float)
    vectors = np.asarray(eigenvectors, dtype=np.complex128)
    if k_grid.ndim != 1 or k_grid.size < 3:
        raise ValueError("k_grid must be one-dimensional with at least three points")
    if vectors.shape[0] != k_grid.size or vectors.ndim != 3 or vectors.shape[1] != vectors.shape[2]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    steps = np.diff(k_grid)
    if not np.allclose(steps, steps[0]):
        raise ValueError("berry_connection currently requires a uniform k grid")

    if smooth_gauge:
        vectors = smooth_eigenvector_phases(vectors)

    dk = float(steps[0])
    derivatives = (np.roll(vectors, -1, axis=0) - np.roll(vectors, 1, axis=0)) / (2.0 * dk)
    xi = np.empty_like(vectors)
    for ik in range(k_grid.size):
        xi[ik] = 1j * vectors[ik].conj().T @ derivatives[ik]
        xi[ik] = 0.5 * (xi[ik] + xi[ik].conj().T)
    return xi


def covariant_derivative_matrix(
    matrix_grid: np.ndarray,
    berry_connection_grid: np.ndarray,
    k_grid: np.ndarray,
) -> np.ndarray:
    """Return [D_k, M] = d_k M - i[xi, M] on a periodic 1D grid."""
    matrices = np.asarray(matrix_grid, dtype=np.complex128)
    xi = np.asarray(berry_connection_grid, dtype=np.complex128)
    k_grid = np.asarray(k_grid, dtype=float)

    if matrices.ndim != 3 or matrices.shape[1] != matrices.shape[2]:
        raise ValueError("matrix_grid must have shape (Nk, Nb, Nb)")
    if xi.shape != matrices.shape:
        raise ValueError("berry_connection_grid must match matrix_grid shape")
    if k_grid.ndim != 1 or k_grid.size != matrices.shape[0]:
        raise ValueError("k_grid size must match matrix_grid")

    steps = np.diff(k_grid)
    if not np.allclose(steps, steps[0]):
        raise ValueError("covariant_derivative_matrix requires a uniform k grid")

    dk = float(steps[0])
    derivative = (np.roll(matrices, -1, axis=0) - np.roll(matrices, 1, axis=0)) / (2.0 * dk)
    commutator = xi @ matrices - matrices @ xi
    return derivative - 1j * commutator
