"""Macroscopic current from RDM trajectories over the first Brillouin zone."""

from __future__ import annotations

from typing import Callable, Mapping

import numpy as np

try:
    from .Band_Solver import build_hamiltonian
    from .Geometry import velocity_matrix
    from .RDM_Common import drifted_k, propagate_single_k_plane_wave_rdm
    from .config import make_k_grid, normalize_params
except ImportError:
    from Band_Solver import build_hamiltonian
    from Geometry import velocity_matrix
    from RDM_Common import drifted_k, propagate_single_k_plane_wave_rdm
    from config import make_k_grid, normalize_params

VelocityFunction = Callable[[float, Mapping], np.ndarray]
RDMSolverFunction = Callable[..., tuple[np.ndarray, np.ndarray]]


def current_for_k_trajectory(
    params: Mapping | None,
    k0: float,
    time_grid: np.ndarray,
    rho_trajectory: np.ndarray,
    velocity_function: VelocityFunction | None = None,
) -> np.ndarray:
    """Return J_k(t) = Re Tr[v(k(t)) rho(t)] for one initial k0."""
    p = normalize_params(params)
    if velocity_function is None:
        velocity_function = velocity_matrix

    contribution = np.empty(time_grid.size, dtype=float)
    for it, t in enumerate(time_grid):
        k_current = drifted_k(float(k0), float(t), p)
        v_current = velocity_function(k_current, p)
        contribution[it] = float(np.trace(v_current @ rho_trajectory[it]).real)
    return contribution


def total_current(
    params: Mapping | None = None,
    k_grid: np.ndarray | None = None,
    velocity_function: VelocityFunction | None = None,
    rdm_solver_function: RDMSolverFunction | None = None,
    k_weight: float | None = None,
    time_span: tuple[float, float] | None = None,
    build_H=build_hamiltonian,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the RDM current over the first Brillouin zone.

    Returns ``(time_grid, total_current)``.  The integration uses the uniform
    weight ``dk/(2*pi)`` for endpoint-excluding grids.
    """
    p = normalize_params(params)
    if k_grid is None:
        k_grid, default_weight = make_k_grid(p, num_k=101)
    else:
        k_grid = np.asarray(k_grid, dtype=float)
        if k_grid.ndim != 1 or k_grid.size == 0:
            raise ValueError("k_grid must be a non-empty one-dimensional array")
        default_weight = p["b"] / k_grid.size / (2.0 * np.pi)

    if k_weight is None:
        k_weight = default_weight
    if velocity_function is None:
        velocity_function = velocity_matrix
    if rdm_solver_function is None:
        rdm_solver_function = propagate_single_k_plane_wave_rdm

    time_grid_ref: np.ndarray | None = None
    accumulated: np.ndarray | None = None

    for k0 in k_grid:
        time_grid, rho_trajectory = rdm_solver_function(
            p,
            float(k0),
            time_span=time_span,
            build_H=build_H,
        )
        if time_grid_ref is None:
            time_grid_ref = time_grid
            accumulated = np.zeros(time_grid.size, dtype=float)
        elif not np.allclose(time_grid_ref, time_grid):
            raise ValueError("RDM solver returned inconsistent time grids")

        accumulated += current_for_k_trajectory(
            p,
            float(k0),
            time_grid,
            rho_trajectory,
            velocity_function=velocity_function,
        ) * float(k_weight)

    assert time_grid_ref is not None and accumulated is not None
    return time_grid_ref, -p["e_charge"] * accumulated


def calculate_current(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for ``total_current``."""
    return total_current(*args, **kwargs)

