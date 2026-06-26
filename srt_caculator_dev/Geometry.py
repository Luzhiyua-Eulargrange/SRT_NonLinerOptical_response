"""Geometry and velocity matrices for the 1D continuum model."""

from __future__ import annotations

from typing import Mapping

import numpy as np

try:
    from .Band_Solver import diagonalize
    from .config import normalize_params
except ImportError:
    from Band_Solver import diagonalize
    from config import normalize_params

def velocity_matrix(
    k_point: float,
    params: Mapping | None = None,
    basis: str = "plane_wave",
) -> np.ndarray:
    """
    Compute the velocity operator matrix at k.
    The plane-wave representation uses v = (1/hbar) dH/dk. If basis="band",
    the matrix is transformed to the instantaneous eigenbasis of H(k).
    """
    p = normalize_params(params)
    ells = p["ell_list"]
    num_g = p["NG"]
    velocity = np.zeros((p["Nb"], p["Nb"]), dtype=np.complex128)

    hbar = p["hbar"]
    prefactor = hbar / p["m"]

    for i, ell in enumerate(ells):
        G = float(ell) * p["b"]
        velocity[i, i] = prefactor * (k_point + G - p["kappa"])
        velocity[num_g + i, num_g + i] = -prefactor * (k_point + G)

    if basis == "plane_wave":
        return velocity
    if basis == "band":
        _, eigenvectors = diagonalize(k_point, p)
        return eigenvectors.conj().T @ velocity @ eigenvectors
    raise ValueError("basis must be 'plane_wave' or 'band'")


def smooth_eigenvector_gauge(eigenvectors: np.ndarray) -> np.ndarray:
    #Align adjacent eigenvector phases along a one-dimensional k grid.
    # Align eigenvector phases along k-grid via parallel transport to enable stable numerical differentiation for Berry connection.
    vectors = np.array(eigenvectors, dtype=np.complex128, copy=True)
    if vectors.ndim != 3 or vectors.shape[1] != vectors.shape[2]:
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")

    for ik in range(1, vectors.shape[0]):
        for band in range(vectors.shape[2]):
            overlap = np.vdot(vectors[ik - 1, :, band], vectors[ik, :, band])
            if abs(overlap) > 0.0:
                vectors[ik, :, band] *= np.exp(-1j * np.angle(overlap))
    return vectors


def berry_connection(
    k_grid: np.ndarray,
    eigenvectors: np.ndarray,
    smooth_gauge: bool = True,
) -> np.ndarray:
    #Compute xi_ss'(k) = i <u_ks | d_k u_ks'> by finite differences.
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
        vectors = smooth_eigenvector_gauge(vectors)

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


def band_velocity_matrices(
    k_grid: np.ndarray,
    params: Mapping | None = None,
    eigenvectors: np.ndarray | None = None,
) -> np.ndarray:
    """
    Return U^dagger (1/hbar dH/dk) U for every k.

    This analytic path is used by the production velocity-gauge evolution. The
    covariant-derivative form from Eq. (28) is used separately as a sparse
    diagnostic because it depends on finite-difference Berry-connection data.
    """
    p = normalize_params(params)
    k_grid = np.asarray(k_grid, dtype=float)
    velocities = np.empty((k_grid.size, p["Nb"], p["Nb"]), dtype=np.complex128)

    if eigenvectors is not None:
        eigenvectors = np.asarray(eigenvectors, dtype=np.complex128)
        if eigenvectors.shape != velocities.shape:
            raise ValueError(f"eigenvectors must have shape {velocities.shape}")

    for ik, k_point in enumerate(k_grid):
        v_plane = velocity_matrix(float(k_point), p, basis="plane_wave")
        if eigenvectors is None:
            velocities[ik] = velocity_matrix(float(k_point), p, basis="band")
        else:
            u = eigenvectors[ik]
            velocities[ik] = u.conj().T @ v_plane @ u
    return velocities
