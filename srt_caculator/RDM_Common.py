"""Shared helpers for reduced-density-matrix gauge solvers."""

from __future__ import annotations

from typing import Callable, Mapping

import numpy as np

try:
    from scipy.integrate import solve_ivp
except ImportError:  # pragma: no cover - used only when scipy is unavailable.
    solve_ivp = None

try:
    from .Band_Solver import build_hamiltonian, diagonalize
    from .config import normalize_params
    from .config import fold_to_fbz
except ImportError:
    from Band_Solver import build_hamiltonian, diagonalize
    from config import normalize_params
    from config import fold_to_fbz

HamiltonianBuilder = Callable[[float, Mapping], np.ndarray]


def smooth_envelope(t: float, params: Mapping | None = None) -> float:
    """Sin-squared pulse envelope with optional smooth switch-on."""
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
    """Vector potential used by the velocity-gauge RDM equation."""
    p = normalize_params(params)
    return float((p["E0"] / p["omega"]) * smooth_envelope(t, p) * np.cos(p["omega"] * t))


def drifted_k(k0: float, t: float, params: Mapping | None = None) -> float:
    """Return k0 + e A(t) / hbar folded to the first Brillouin zone."""
    p = normalize_params(params)
    k_drift = float(k0) + p["e_charge"] * vector_potential(t, p) / p["hbar"]
    return float(fold_to_fbz(k_drift, p["b"]))


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


def initial_plane_wave_density_matrix(
    k0: float,
    params: Mapping | None = None,
    build_H: HamiltonianBuilder | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct the zero-temperature adiabatic RDM in the plane-wave basis."""
    p = normalize_params(params)
    if build_H is None:
        energies, eigenvectors = diagonalize(float(k0), p)
    else:
        hamiltonian = build_H(float(k0), p)
        energies, eigenvectors = np.linalg.eigh(hamiltonian)

    occupations = (energies < p["fermi_energy"]).astype(float)
    rho0 = eigenvectors @ np.diag(occupations) @ eigenvectors.conj().T
    return rho0.astype(np.complex128), energies, eigenvectors


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


def propagate_single_k_plane_wave_rdm(
    params: Mapping | None,
    k0: float,
    time_span: tuple[float, float] | None = None,
    build_H: HamiltonianBuilder | None = None,
    time_grid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate the legacy single-k RDM in the plane-wave basis.

    This keeps the old current-integration path available from the shared
    module. New gauge-resolved calculations should use RDM_Velocity_Gauge or
    RDM_Length_Gauge instead.
    """
    p = normalize_params(params)
    if build_H is None:
        build_H = build_hamiltonian

    t_grid = make_time_grid(p, time_span=time_span, time_grid=time_grid)
    rho0, _, _ = initial_plane_wave_density_matrix(float(k0), p, build_H=build_H)
    shape = rho0.shape

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        rho = y.reshape(shape)
        k_current = drifted_k(float(k0), t, p)
        hamiltonian = build_H(k_current, p)
        commutator = hamiltonian @ rho - rho @ hamiltonian
        return (-1j * commutator / p["hbar"]).reshape(-1)

    values = solve_complex_ode(
        rhs,
        rho0.reshape(-1),
        t_grid,
        p,
        error_label="single-k plane-wave RDM propagation",
    )
    rho_trajectory = values.reshape((t_grid.size,) + shape)
    return t_grid, hermitize_density_trajectory(rho_trajectory)


def propagate_rdm(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for propagate_single_k_plane_wave_rdm."""
    return propagate_single_k_plane_wave_rdm(*args, **kwargs)


def solve_rdm(*args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for propagate_single_k_plane_wave_rdm."""
    return propagate_single_k_plane_wave_rdm(*args, **kwargs)
