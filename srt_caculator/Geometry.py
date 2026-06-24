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
    """Compute the velocity operator matrix at k.

    The plane-wave representation uses ``v = (1/hbar) dH/dk``.  If
    ``basis='band'``, the matrix is transformed to the instantaneous eigenbasis
    of ``H(k)``.
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


def berry_connection(*_args, **_kwargs) -> np.ndarray:
    """Placeholder for future covariant-derivative geometry.

    The current RDM solver only requires the velocity matrix.  This function is
    intentionally explicit so callers do not silently get an incorrect Berry
    connection approximation.
    """
    raise NotImplementedError("Berry connection is not implemented in the base solver")
