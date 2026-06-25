"""Length-gauge RDM propagation based on Phys. Rev. B 96, 035431 Eq. (55)."""

from __future__ import annotations

from typing import Mapping

import numpy as np

try:
    from .Geometry import berry_connection, covariant_derivative_matrix
    from .RDM_Common import (
        electric_field,
        hermitize_density_trajectory,
        initial_band_density_matrices,
        make_time_grid,
        solve_complex_ode,
    )
    from .config import normalize_params
except ImportError:
    from Geometry import berry_connection, covariant_derivative_matrix
    from RDM_Common import (
        electric_field,
        hermitize_density_trajectory,
        initial_band_density_matrices,
        make_time_grid,
        solve_complex_ode,
    )
    from config import normalize_params


def _validate_band_inputs(
    k_grid: np.ndarray,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k_grid = np.asarray(k_grid, dtype=float)
    energies = np.asarray(energies, dtype=float)
    eigenvectors = np.asarray(eigenvectors, dtype=np.complex128)
    if k_grid.ndim != 1:
        raise ValueError("k_grid must be one-dimensional")
    if energies.ndim != 2 or energies.shape[0] != k_grid.size:
        raise ValueError("energies must have shape (Nk, Nb)")
    if eigenvectors.shape != (k_grid.size, energies.shape[1], energies.shape[1]):
        raise ValueError("eigenvectors must have shape (Nk, Nb, Nb)")
    return k_grid, energies, eigenvectors


def length_gauge_rhs(
    t: float,
    rho_grid: np.ndarray,
    params: Mapping | None,
    k_grid: np.ndarray,
    energies: np.ndarray,
    berry_connection_grid: np.ndarray,
) -> np.ndarray:
    """Return d rho / dt from Eq. (55) on a k grid.

    Eq. (55): (i hbar d_t - epsilon_ss') rho = i e E(t) [D, rho].
    In 1D, [D, rho] = d_k rho - i [xi, rho].
    """
    p = normalize_params(params)
    rho_grid = np.asarray(rho_grid, dtype=np.complex128)
    energy_diff = energies[:, :, None] - energies[:, None, :]
    covariant_derivative = covariant_derivative_matrix(rho_grid, berry_connection_grid, k_grid)
    coherent_term = -1j * energy_diff * rho_grid / p["hbar"]
    field_term = p["e_charge"] * electric_field(t, p) * covariant_derivative / p["hbar"]
    return coherent_term + field_term


def propagate_length_gauge_rdm(
    params: Mapping | None,
    k_grid: np.ndarray,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
    time_span: tuple[float, float] | None = None,
    time_grid: np.ndarray | None = None,
    berry_connection_grid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate length-gauge RDMs for all k points.

    Returns (time_grid, rho_trajectory), where rho_trajectory has shape
    (Nt, Nk, Nb, Nb) in the band basis.
    """
    p = normalize_params(params)
    k_grid, energies, eigenvectors = _validate_band_inputs(k_grid, energies, eigenvectors)
    if berry_connection_grid is None:
        berry_connection_grid = berry_connection(k_grid, eigenvectors)
    else:
        berry_connection_grid = np.asarray(berry_connection_grid, dtype=np.complex128)

    if berry_connection_grid.shape != (k_grid.size, energies.shape[1], energies.shape[1]):
        raise ValueError("berry_connection_grid must have shape (Nk, Nb, Nb)")

    t_grid = make_time_grid(p, time_span=time_span, time_grid=time_grid)
    rho0 = initial_band_density_matrices(energies, p)
    shape = rho0.shape

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        rho = y.reshape(shape)
        return length_gauge_rhs(t, rho, p, k_grid, energies, berry_connection_grid).reshape(-1)

    values = solve_complex_ode(
        rhs,
        rho0.reshape(-1),
        t_grid,
        p,
        error_label="length-gauge RDM propagation",
    )
    rho_trajectory = values.reshape((t_grid.size,) + shape)
    return t_grid, hermitize_density_trajectory(rho_trajectory)


def solve_rdm_length_gauge(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for propagate_length_gauge_rdm."""
    return propagate_length_gauge_rdm(*args, **kwargs)
