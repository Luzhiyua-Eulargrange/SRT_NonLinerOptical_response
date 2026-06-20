"""Non-perturbative reduced density matrix propagation."""

from __future__ import annotations

from typing import Callable, Mapping

import numpy as np

try:
    from scipy.integrate import solve_ivp
except ImportError:  # pragma: no cover - used only when scipy is unavailable.
    solve_ivp = None

try:
    from .Band_Solver import build_hamiltonian, diagonalize
    from .config import fold_to_fbz, normalize_params
except ImportError:
    from Band_Solver import build_hamiltonian, diagonalize
    from config import fold_to_fbz, normalize_params

HamiltonianBuilder = Callable[[float, Mapping], np.ndarray]


def smooth_envelope(t: float, params: Mapping | None = None) -> float:
    """Sin-squared pulse envelope with optional switch-on time."""
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


def vector_potential(t: float, params: Mapping | None = None) -> float:
    """Vector potential compatible with the module requirements.

    This closed-form approximation follows the requested convention
    ``A(t) = E0 / omega * envelope(t) * cos(omega t)``.  For slowly varying
    envelopes this gives ``E(t) ~= -dA/dt`` up to envelope-derivative terms.
    """
    p = normalize_params(params)
    return float((p["E0"] / p["omega"]) * smooth_envelope(t, p) * np.cos(p["omega"] * t))


def drifted_k(k0: float, t: float, params: Mapping | None = None) -> float:
    """Return k0 + e A(t) / hbar folded to the first Brillouin zone."""
    p = normalize_params(params)
    k_drift = float(k0) + p["e_charge"] * vector_potential(t, p) / p["hbar"]
    return float(fold_to_fbz(k_drift, p["b"]))


def initial_density_matrix(
    k0: float,
    params: Mapping | None = None,
    build_H: HamiltonianBuilder | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct the zero-temperature adiabatic ground-state RDM."""
    p = normalize_params(params)

    if build_H is None:
        energies, eigenvectors = diagonalize(k0, p)
    else:
        H0 = build_H(float(k0), p)
        energies, eigenvectors = np.linalg.eigh(H0)

    occupations = (energies < p["fermi_energy"]).astype(float)
    rho0 = eigenvectors @ np.diag(occupations) @ eigenvectors.conj().T
    return rho0.astype(np.complex128), energies, eigenvectors


def _rhs(
    t: float,
    y: np.ndarray,
    k0: float,
    params: Mapping,
    build_H: HamiltonianBuilder,
) -> np.ndarray:
    p = normalize_params(params)
    nb = p["Nb"]
    rho = y.reshape((nb, nb))
    k_current = drifted_k(k0, t, p)
    H_current = build_H(k_current, p)
    commutator = H_current @ rho - rho @ H_current
    drho_dt = -1j * commutator / p["hbar"]
    return drho_dt.reshape(-1)


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


def propagate_rdm(
    params: Mapping | None,
    k0: float,
    time_span: tuple[float, float] | None = None,
    build_H: HamiltonianBuilder | None = None,
    time_grid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate the RDM for one initial momentum k0.

    Returns ``(time_grid, rho_trajectory)`` where rho has shape ``(Nt, Nb, Nb)``.
    """
    p = normalize_params(params)
    if build_H is None:
        build_H = build_hamiltonian

    if time_span is None:
        time_span = (0.0, p["pulse_duration"])
    t_start, t_end = float(time_span[0]), float(time_span[1])
    if t_end <= t_start:
        raise ValueError("time_span must satisfy end > start")

    if time_grid is None:
        time_grid = np.linspace(t_start, t_end, p["num_time_steps"])
    else:
        time_grid = np.asarray(time_grid, dtype=float)
        if time_grid.ndim != 1 or time_grid.size < 2:
            raise ValueError("time_grid must be a one-dimensional array with at least two points")

    rho0, _, _ = initial_density_matrix(k0, p, build_H)
    y0 = rho0.reshape(-1)

    rhs = lambda t, y: _rhs(t, y, float(k0), p, build_H)
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
            raise RuntimeError(f"RDM propagation failed: {solution.message}")
        values = solution.y.T
    else:
        values = _rk4_fallback(rhs, y0, time_grid)

    rho_trajectory = values.reshape((time_grid.size, p["Nb"], p["Nb"]))
    rho_trajectory = 0.5 * (rho_trajectory + rho_trajectory.conj().transpose(0, 2, 1))
    return time_grid, rho_trajectory


def solve_rdm(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for ``propagate_rdm``."""
    return propagate_rdm(*args, **kwargs)

