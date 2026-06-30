"""Debug utilities for inspecting solver parameters and runtime state."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
from typing import Mapping

import numpy as np

from config import DEFAULT_PARAMS


def print_params(params: Mapping, title: str = "PARAMS") -> None:
    """Print a parameter mapping in a readable key-value format."""
    print(f"{title}:")
    for key, value in params.items():
        print(f"  {key}: {value}")


def print_default_params() -> None:
    """Print the default parameter list used by the solver."""
    print_params(DEFAULT_PARAMS, title="DEFAULT_PARAMS")


def plot_bands_and_print_eigenvectors(
    k_grid: np.ndarray,
    energies: np.ndarray,
    eigenvectors: np.ndarray,
    output_path: str = "band_structure.png",
    k_indices: Iterable[int] = (0,),
    print_eigenvectors: bool = True,
) -> None:
    """Plot band energies and print selected eigenvector matrices."""
    k_grid = np.asarray(k_grid, dtype=float)
    energies = np.asarray(energies, dtype=float)
    eigenvectors = np.asarray(eigenvectors, dtype=np.complex128)

    if k_grid.ndim != 1:
        raise ValueError("k_grid must be one-dimensional")
    if energies.ndim != 2:
        raise ValueError("energies must be two-dimensional with shape (Nk, Nb)")
    if eigenvectors.ndim != 3:
        raise ValueError("eigenvectors must be three-dimensional with shape (Nk, Nb, Nb)")
    if energies.shape[0] != k_grid.size:
        raise ValueError("energies and k_grid have inconsistent Nk")
    expected_shape = (k_grid.size, energies.shape[1], energies.shape[1])
    if eigenvectors.shape != expected_shape:
        raise ValueError(f"eigenvectors shape must be {expected_shape}, got {eigenvectors.shape}")

    print("Band result debug check:")
    print("  k_grid shape:", k_grid.shape)
    print("  energies shape:", energies.shape)
    print("  eigenvectors shape:", eigenvectors.shape)
    print("  energies finite:", bool(np.all(np.isfinite(energies))))

    identity = np.eye(energies.shape[1], dtype=np.complex128)
    max_orthogonality_error = 0.0
    for matrix in eigenvectors:
        error = np.max(np.abs(matrix.conj().T @ matrix - identity))
        max_orthogonality_error = max(max_orthogonality_error, float(error))
    print("  max eigenvector orthogonality error:", max_orthogonality_error)

    if print_eigenvectors:
        with np.printoptions(precision=6, suppress=True):
            for index in k_indices:
                if index < 0:
                    index += k_grid.size
                if index < 0 or index >= k_grid.size:
                    raise IndexError(f"k index {index} is out of range for Nk={k_grid.size}")
                print(f"Eigenvectors at k_grid[{index}] = {k_grid[index]}:")
                print(eigenvectors[index])

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped band plot.")
        return

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for band_index in range(energies.shape[1]):
        ax.plot(k_grid, energies[:, band_index], linewidth=1.2)
    ax.set_xlabel("k")
    ax.set_ylabel("Energy")
    ax.set_title("Band structure")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved band plot to {output.resolve()}")


def debug_rdm_trajectory(
    time_grid: np.ndarray,
    rho_trajectory: np.ndarray,
    title: str = "RDM trajectory",
) -> None:
    """Print basic consistency checks for an RDM trajectory."""
    time_grid = np.asarray(time_grid, dtype=float)
    rho = np.asarray(rho_trajectory, dtype=np.complex128)
    if time_grid.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if rho.shape[0] != time_grid.size:
        raise ValueError("rho_trajectory first axis must match time_grid")
    if rho.ndim not in (3, 4) or rho.shape[-1] != rho.shape[-2]:
        raise ValueError("rho_trajectory must have shape (Nt, Nb, Nb) or (Nt, Nk, Nb, Nb)")

    hermitian_error = np.max(np.abs(rho - rho.conj().swapaxes(-1, -2)))
    # For closed-system commutator dynamics, Tr(rho) is conserved and tracks
    # the total occupation. A large drift usually means an equation or ODE issue.
    traces = np.trace(rho, axis1=-2, axis2=-1)
    trace_drift = np.max(np.abs(traces - traces[0]))

    print(f"{title}:")
    print("  time_grid shape:", time_grid.shape)
    print("  rho_trajectory shape:", rho.shape)
    print("  max Hermitian error:", float(hermitian_error))
    print("  max trace drift:", float(trace_drift))
