"""Velocity-gauge RDM propagation based on Phys. Rev. B 96, 035431 Eq. (56)."""

from __future__ import annotations

from typing import Mapping

import numpy as np

try:
    from .Geometry import band_velocity_matrices
    from .RDM_Common import (
        hermitize_density_trajectory,
        initial_band_density_matrices,
        make_time_grid,
        solve_complex_ode,
        vector_potential,
    )
    from .config import normalize_params
except ImportError:
    from Geometry import band_velocity_matrices
    from RDM_Common import (
        hermitize_density_trajectory,
        initial_band_density_matrices,
        make_time_grid,
        solve_complex_ode,
        vector_potential,
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


def velocity_gauge_rhs(
    t: float,
    rho_grid: np.ndarray,
    params: Mapping | None,
    energies: np.ndarray,
    velocity_grid: np.ndarray,
) -> np.ndarray:
    """Return d rho / dt from Eq. (56) on a k grid.

    Eq. (56): (i hbar d_t - epsilon_ss') rho = e A(t) [v, rho].
    """
    p = normalize_params(params)
    rho_grid = np.asarray(rho_grid, dtype=np.complex128)
    energy_diff = energies[:, :, None] - energies[:, None, :]
    commutator = velocity_grid @ rho_grid - rho_grid @ velocity_grid
    rhs_operator = energy_diff * rho_grid + p["e_charge"] * vector_potential(t, p) * commutator
    return -1j * rhs_operator / p["hbar"]


def propagate_velocity_gauge_rdm(
    params: Mapping | None,
    k_grid: np.ndarray,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
    time_span: tuple[float, float] | None = None,
    time_grid: np.ndarray | None = None,
    velocity_grid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate velocity-gauge RDMs for all k points.

    Returns (time_grid, rho_trajectory), where rho_trajectory has shape
    (Nt, Nk, Nb, Nb) in the band basis.
    """
    p = normalize_params(params)
    k_grid, energies, eigenvectors = _validate_band_inputs(k_grid, energies, eigenvectors)
    if velocity_grid is None:
        velocity_grid = band_velocity_matrices(k_grid, p, eigenvectors=eigenvectors)
    else:
        velocity_grid = np.asarray(velocity_grid, dtype=np.complex128)

    if velocity_grid.shape != (k_grid.size, energies.shape[1], energies.shape[1]):
        raise ValueError("velocity_grid must have shape (Nk, Nb, Nb)")

    t_grid = make_time_grid(p, time_span=time_span, time_grid=time_grid)
    rho0 = initial_band_density_matrices(energies, p)
    shape = rho0.shape

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        rho = y.reshape(shape)
        return velocity_gauge_rhs(t, rho, p, energies, velocity_grid).reshape(-1)

    values = solve_complex_ode(
        rhs,
        rho0.reshape(-1),
        t_grid,
        p,
        error_label="velocity-gauge RDM propagation",
    )
    rho_trajectory = values.reshape((t_grid.size,) + shape)
    return t_grid, hermitize_density_trajectory(rho_trajectory)


def solve_rdm_velocity_gauge(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for propagate_velocity_gauge_rdm."""
    return propagate_velocity_gauge_rdm(*args, **kwargs)
