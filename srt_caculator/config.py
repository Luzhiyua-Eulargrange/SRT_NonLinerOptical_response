"""
Configuration helpers for the SRT RDM solver.

The project uses a dictionary application programming interface so scripts can override only the parameters they need.
'normalize_params' fills defaults, validates inputs, and adds derived values such as the reciprocal lattice vector and matrix size.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class DefaultParams:
    # Static one-dimensional two-component continuum model.
    a0: float = 1.0
    m: float = 1.0
    hbar: float = 1.0
    kappa: float = 0.15
    v0: float = 0.2
    vA: complex = 0.08 + 0.02j
    vB: complex =-0.06 + 0.01j
    w1: complex = 0.15 + 0.00j
    w2: complex = 0.10 - 0.03j
    L: int = 3  #Reciprocal lattice vector truncation parameter

    # Field and time propagation parameters.
    e_charge: float = 1.0
    E0: float = 0.02
    omega: float = 0.5
    pulse_duration: float = 40.0
    t_switch: float = 8.0
    num_time_steps: int = 401
    fermi_energy: float = 0.0
    rtol: float = 1e-8  #relative tolerance
    atol: float = 1e-10 #absolute tolerance


DEFAULT_PARAMS = asdict(DefaultParams())


def normalize_params(params: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Return a validated parameter dictionary with derived quantities added."""
    out: dict[str, Any] = dict(DEFAULT_PARAMS)
    if params is not None:
        out.update(dict(params))

    out["a0"] = float(out["a0"])
    out["m"] = float(out["m"])
    out["hbar"] = float(out["hbar"])
    out["kappa"] = float(out["kappa"])
    out["v0"] = float(out["v0"])
    out["vA"] = complex(out["vA"])
    out["vB"] = complex(out["vB"])
    out["w1"] = complex(out["w1"])
    out["w2"] = complex(out["w2"])
    out["L"] = int(out["L"])
    out["e_charge"] = float(out["e_charge"])
    out["E0"] = float(out["E0"])
    out["omega"] = float(out["omega"])
    out["pulse_duration"] = float(out["pulse_duration"])
    out["t_switch"] = float(out["t_switch"])
    out["num_time_steps"] = int(out["num_time_steps"])
    out["fermi_energy"] = float(out["fermi_energy"])
    out["rtol"] = float(out["rtol"])
    out["atol"] = float(out["atol"])

    if out["a0"] <= 0.0:
        raise ValueError("a0 must be positive")
    if out["m"] <= 0.0:
        raise ValueError("m must be positive")
    if out["hbar"] <= 0.0:
        raise ValueError("hbar must be positive")
    if out["L"] < 0:
        raise ValueError("L must be non-negative")
    if out["omega"] == 0.0:
        raise ValueError("omega must be non-zero")
    if out["pulse_duration"] <= 0.0:
        raise ValueError("pulse_duration must be positive")
    if out["t_switch"] < 0.0:
        raise ValueError("t_switch must be non-negative")
    if out["num_time_steps"] < 2:
        raise ValueError("num_time_steps must be at least 2")

    out["b"] = 2.0 * np.pi / out["a0"]
    out["ell_list"] = np.arange(-out["L"], out["L"] + 1, dtype=int)
    out["NG"] = int(2 * out["L"] + 1)
    out["Nb"] = int(2 * out["NG"])
    return out


def fold_to_fbz(k_value: float | np.ndarray, b: float) -> float | np.ndarray:
    """Fold momentum to the interval [-b/2, b/2)."""
    return (np.asarray(k_value) + 0.5 * b) % b - 0.5 * b


def make_k_grid(params: Mapping[str, Any] | None = None, num_k: int = 101) -> tuple[np.ndarray, float]:
    """Return an endpoint-excluding FBZ grid and integration weight dk/(2*pi) for num integration."""
    p = normalize_params(params)
    if num_k < 2:
        raise ValueError("num_k must be at least 2")
    k_grid = np.linspace(-0.5 * p["b"], 0.5 * p["b"], int(num_k), endpoint=False)
    dk = p["b"] / int(num_k)
    return k_grid, dk / (2.0 * np.pi)
