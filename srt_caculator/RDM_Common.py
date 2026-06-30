"""Shared helpers for reduced-density-matrix solvers."""

from __future__ import annotations

from typing import Callable, Mapping

import numpy as np

try:
    from scipy.integrate import solve_ivp
except ImportError:  # pragma: no cover - used only when scipy is unavailable.
    solve_ivp = None

try:
    from .config import normalize_params
except ImportError:
    from config import normalize_params


def smooth_envelope(t: float, params: Mapping | None = None) -> float:
    # Sin-squared pulse envelope with optional smooth switch-on.
    p = normalize_params(params)
    duration = p["pulse_duration"]
    if t < 0.0 or t > duration:
        return 0.0

    base = np.sin(np.pi * t / duration) ** 2
    if p["t_switch"] <= 0.0:
        return float(base)

    ramp = min(1.0, max(0.0, t / p["t_switch"]))
    return float(base * ramp * ramp * (3.0 - 2.0 * ramp))


def electric_field(t: float, params: Mapping | None = None) -> float:
    """Electric field E(t) = E0 envelope(t) sin(omega t)."""
    p = normalize_params(params)
    return float(p["E0"] * smooth_envelope(t, p) * np.sin(p["omega"] * t))


def make_time_grid(
    params: Mapping | None = None,
    time_span: tuple[float, float] | None = None,
    time_grid: np.ndarray | None = None,
) -> np.ndarray:
    """Return a validated time grid."""
    p = normalize_params(params)
    if time_grid is not None:
        out = np.asarray(time_grid, dtype=float)
        if out.ndim != 1 or out.size < 2:
            raise ValueError("time_grid must be a one-dimensional array with at least two points")
        return out

    if time_span is None:
        time_span = (0.0, p["pulse_duration"])
    t_start, t_end = float(time_span[0]), float(time_span[1])
    if t_end <= t_start:
        raise ValueError("time_span must satisfy end > start")
    return np.linspace(t_start, t_end, p["num_time_steps"])


def initial_band_density_matrices(
    energies: np.ndarray,
    params: Mapping | None = None,
) -> np.ndarray:
    """Build zero-temperature diagonal RDMs in the band basis."""
    p = normalize_params(params)
    energies = np.asarray(energies, dtype=float)
    if energies.ndim != 2:
        raise ValueError("energies must have shape (Nk, Nb)")

    occupations = (energies < p["fermi_energy"]).astype(float)
    nb = energies.shape[1]
    rho = np.zeros((energies.shape[0], nb, nb), dtype=np.complex128)
    diag = np.arange(nb)
    rho[:, diag, diag] = occupations
    return rho


def hermitize_density_trajectory(rho_trajectory: np.ndarray) -> np.ndarray:
    """Remove small anti-Hermitian numerical drift from RDM trajectories."""
    return 0.5 * (rho_trajectory + rho_trajectory.conj().swapaxes(-1, -2))


def _rk4_fallback(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    time_grid: np.ndarray,
) -> np.ndarray:
    values = np.empty((time_grid.size, y0.size), dtype=np.complex128)
    values[0] = y0
    for i in range(time_grid.size - 1):
        t = float(time_grid[i])
        dt = float(time_grid[i + 1] - time_grid[i])
        y = values[i]
        k1 = rhs(t, y)
        k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
        k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
        k4 = rhs(t + dt, y + dt * k3)
        values[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return values


def solve_complex_ode(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    time_grid: np.ndarray,
    params: Mapping | None = None,
    error_label: str = "RDM propagation",
) -> np.ndarray:
    """Integrate a complex-valued ODE and return values on time_grid."""
    p = normalize_params(params)
    y0 = np.asarray(y0, dtype=np.complex128)
    time_grid = np.asarray(time_grid, dtype=float)

    if solve_ivp is not None:
        solution = solve_ivp(
            rhs,
            (float(time_grid[0]), float(time_grid[-1])),
            y0,
            method="RK45",
            t_eval=time_grid,
            rtol=p["rtol"],
            atol=p["atol"],
        )
        if not solution.success:
            raise RuntimeError(f"{error_label} failed: {solution.message}")
        return solution.y.T

    return _rk4_fallback(rhs, y0, time_grid)
