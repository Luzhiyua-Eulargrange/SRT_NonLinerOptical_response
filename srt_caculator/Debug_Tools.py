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


def load_npz_result(filename: str | Path) -> dict:
    """Load an npz result file into a plain dictionary."""
    data = np.load(filename)
    return {key: data[key] for key in data.files}


def nearest_grid_indices(grid: np.ndarray, sample_points: Iterable[float]) -> np.ndarray:
    """Return indices of grid points closest to the requested sample points."""
    values = np.asarray(grid, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("grid must be a non-empty one-dimensional array")

    indices = []
    for point in sample_points:
        indices.append(int(np.argmin(np.abs(values - float(point)))))
    return np.asarray(indices, dtype=int)


def default_sample_indices(grid_size: int, num_samples: int = 5) -> np.ndarray:
    """Return stable, unique sample indices across a one-dimensional grid."""
    if int(grid_size) < 1:
        raise ValueError("grid_size must be positive")
    count = max(1, min(int(num_samples), int(grid_size)))
    return np.unique(np.linspace(0, int(grid_size) - 1, count, dtype=int))


def plot_density_time_traces(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    density_total: np.ndarray,
    density_components: np.ndarray | None = None,
    output_path: str | Path = "density_time_traces.png",
    x_indices: Iterable[int] | None = None,
    x_points: Iterable[float] | None = None,
    title: str = "Real-space density time traces",
    subtract_initial: bool = False,
) -> None:
    """Plot density time traces at selected real-space grid points."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(density_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("density_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    components = None
    if density_components is not None:
        components = np.asarray(density_components, dtype=float)
        if components.shape != (time.size, 2, x_values.size):
            raise ValueError("density_components must have shape (Nt, 2, Nx)")
        if subtract_initial:
            components = components - components[0][None, :, :]

    if x_points is not None:
        selected = nearest_grid_indices(x_values, x_points)
    elif x_indices is not None:
        selected = np.asarray(list(x_indices), dtype=int)
    else:
        selected = default_sample_indices(x_values.size, num_samples=5)

    if selected.ndim != 1 or selected.size == 0:
        raise ValueError("at least one x index is required")
    if np.any(selected < 0) or np.any(selected >= x_values.size):
        raise IndexError("x index out of range")
    selected = np.unique(selected)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped density time-trace plot.")
        return

    if components is None:
        fig, axes = plt.subplots(1, 1, figsize=(8.0, 4.8), squeeze=False)
        axes_to_plot = [(axes[0, 0], total, "total density")]
    else:
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True)
        axes_to_plot = [
            (axes[0], total, "total density"),
            (axes[1], components[:, 0, :], "A density"),
            (axes[2], components[:, 1, :], "B density"),
        ]

    for ax, values, ylabel in axes_to_plot:
        for index in selected:
            ax.plot(time, values[:, index], linewidth=1.4, label=f"x={x_values[index]:.4g}")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    axes_to_plot[-1][0].set_xlabel("time")
    fig.suptitle(title)
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved density time-trace plot to {output.resolve()}")


def plot_density_spacetime_map(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    density_total: np.ndarray,
    output_path: str | Path = "density_spacetime_total.png",
    title: str = "Total density n(x,t)",
    subtract_initial: bool = False,
) -> None:
    """Plot a two-dimensional density map with x on the horizontal axis and time on the vertical axis."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(density_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("density_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped density spacetime plot.")
        return

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    if subtract_initial:
        vmax = float(np.max(np.abs(total)))
        color_limit = vmax if vmax > 0.0 else 1.0
        mesh = ax.pcolormesh(
            x_values,
            time,
            total,
            shading="auto",
            cmap="RdBu_r",
            vmin=-color_limit,
            vmax=color_limit,
        )
    else:
        mesh = ax.pcolormesh(x_values, time, total, shading="auto", cmap="viridis")
    ax.set_xlabel("x")
    ax.set_ylabel("time")
    ax.set_title(title)
    colorbar = fig.colorbar(mesh, ax=ax)
    colorbar.set_label("density")
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved density spacetime plot to {output.resolve()}")


def plot_density_time_slices(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    density_total: np.ndarray,
    output_path: str | Path = "density_time_slices.png",
    time_indices: Iterable[int] | None = None,
    time_points: Iterable[float] | None = None,
    title: str = "Density profiles at selected times",
    subtract_initial: bool = False,
) -> None:
    """Plot density profiles n(x,t_i) at selected time points."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(density_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("density_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    if time_points is not None:
        selected = nearest_grid_indices(time, time_points)
    elif time_indices is not None:
        selected = np.asarray(list(time_indices), dtype=int)
    else:
        selected = default_sample_indices(time.size, num_samples=5)

    if selected.ndim != 1 or selected.size == 0:
        raise ValueError("at least one time index is required")
    if np.any(selected < 0) or np.any(selected >= time.size):
        raise IndexError("time index out of range")
    selected = np.unique(selected)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped density time-slice plot.")
        return

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for index in selected:
        ax.plot(x_values, total[index], linewidth=1.5, label=f"t={time[index]:.4g}")
    ax.set_xlabel("x")
    ax.set_ylabel("total density")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved density time-slice plot to {output.resolve()}")


def plot_density_time_traces_from_npz(
    input_path: str | Path = "length_density_result.npz",
    output_path: str | Path = "density_time_traces.png",
    x_indices: Iterable[int] | None = None,
    x_points: Iterable[float] | None = None,
    subtract_initial: bool = False,
) -> None:
    """Load a density npz result and plot time traces at selected x points."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "density_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"density result is missing keys: {missing}")

    plot_density_time_traces(
        result["time_grid"],
        result["x_grid"],
        result["density_total"],
        density_components=result.get("density_components"),
        output_path=output_path,
        x_indices=x_indices,
        x_points=x_points,
        title=f"Density time traces from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_density_spacetime_map_from_npz(
    input_path: str | Path = "length_density_result.npz",
    output_path: str | Path = "density_spacetime_total.png",
    subtract_initial: bool = False,
) -> None:
    """Load a density npz result and plot the total-density spacetime map."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "density_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"density result is missing keys: {missing}")

    plot_density_spacetime_map(
        result["time_grid"],
        result["x_grid"],
        result["density_total"],
        output_path=output_path,
        title=f"Total density from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_density_time_slices_from_npz(
    input_path: str | Path = "length_density_result.npz",
    output_path: str | Path = "density_time_slices.png",
    time_indices: Iterable[int] | None = None,
    time_points: Iterable[float] | None = None,
    subtract_initial: bool = False,
) -> None:
    """Load a density npz result and plot total-density profiles at selected times."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "density_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"density result is missing keys: {missing}")

    plot_density_time_slices(
        result["time_grid"],
        result["x_grid"],
        result["density_total"],
        output_path=output_path,
        time_indices=time_indices,
        time_points=time_points,
        title=f"Density profiles from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_density_standard_views_from_npz(
    input_path: str | Path = "length_density_result.npz",
    output_prefix: str | Path = "density",
    subtract_initial: bool = False,
) -> None:
    """Create the standard static density views: spacetime map, time slices, and time traces."""
    prefix = Path(output_prefix)
    plot_density_spacetime_map_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_spacetime_total.png"),
        subtract_initial=subtract_initial,
    )
    plot_density_time_slices_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_time_slices.png"),
        subtract_initial=subtract_initial,
    )
    plot_density_time_traces_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_time_traces.png"),
        subtract_initial=subtract_initial,
    )


def plot_current_spacetime_map(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    current_total: np.ndarray,
    output_path: str | Path = "current_spacetime_total.png",
    title: str = "Total current j(x,t)",
    subtract_initial: bool = False,
) -> None:
    """Plot a two-dimensional current map with x on the horizontal axis and time on the vertical axis."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(current_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("current_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped current spacetime plot.")
        return

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    vmax = float(np.max(np.abs(total)))
    color_limit = vmax if vmax > 0.0 else 1.0
    mesh = ax.pcolormesh(
        x_values,
        time,
        total,
        shading="auto",
        cmap="RdBu_r",
        vmin=-color_limit,
        vmax=color_limit,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("time")
    ax.set_title(title)
    colorbar = fig.colorbar(mesh, ax=ax)
    colorbar.set_label("current")
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved current spacetime plot to {output.resolve()}")


def plot_current_time_slices(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    current_total: np.ndarray,
    output_path: str | Path = "current_time_slices.png",
    time_indices: Iterable[int] | None = None,
    time_points: Iterable[float] | None = None,
    title: str = "Current profiles at selected times",
    subtract_initial: bool = False,
) -> None:
    """Plot current profiles j(x,t_i) at selected time points."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(current_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("current_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    if time_points is not None:
        selected = nearest_grid_indices(time, time_points)
    elif time_indices is not None:
        selected = np.asarray(list(time_indices), dtype=int)
    else:
        selected = default_sample_indices(time.size, num_samples=5)

    if selected.ndim != 1 or selected.size == 0:
        raise ValueError("at least one time index is required")
    if np.any(selected < 0) or np.any(selected >= time.size):
        raise IndexError("time index out of range")
    selected = np.unique(selected)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped current time-slice plot.")
        return

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for index in selected:
        ax.plot(x_values, total[index], linewidth=1.5, label=f"t={time[index]:.4g}")
    ax.set_xlabel("x")
    ax.set_ylabel("total current")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved current time-slice plot to {output.resolve()}")


def plot_current_time_traces(
    time_grid: np.ndarray,
    x_grid: np.ndarray,
    current_total: np.ndarray,
    current_components: np.ndarray | None = None,
    output_path: str | Path = "current_time_traces.png",
    x_indices: Iterable[int] | None = None,
    x_points: Iterable[float] | None = None,
    title: str = "Real-space current time traces",
    subtract_initial: bool = False,
) -> None:
    """Plot current time traces at selected real-space grid points."""
    time = np.asarray(time_grid, dtype=float)
    x_values = np.asarray(x_grid, dtype=float)
    total = np.asarray(current_total, dtype=float)

    if time.ndim != 1:
        raise ValueError("time_grid must be one-dimensional")
    if x_values.ndim != 1:
        raise ValueError("x_grid must be one-dimensional")
    if total.shape != (time.size, x_values.size):
        raise ValueError("current_total must have shape (Nt, Nx)")
    if subtract_initial:
        total = total - total[0]

    components = None
    if current_components is not None:
        components = np.asarray(current_components, dtype=float)
        if components.shape != (time.size, 2, x_values.size):
            raise ValueError("current_components must have shape (Nt, 2, Nx)")
        if subtract_initial:
            components = components - components[0][None, :, :]

    if x_points is not None:
        selected = nearest_grid_indices(x_values, x_points)
    elif x_indices is not None:
        selected = np.asarray(list(x_indices), dtype=int)
    else:
        selected = default_sample_indices(x_values.size, num_samples=5)

    if selected.ndim != 1 or selected.size == 0:
        raise ValueError("at least one x index is required")
    if np.any(selected < 0) or np.any(selected >= x_values.size):
        raise IndexError("x index out of range")
    selected = np.unique(selected)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped current time-trace plot.")
        return

    if components is None:
        fig, axes = plt.subplots(1, 1, figsize=(8.0, 4.8), squeeze=False)
        axes_to_plot = [(axes[0, 0], total, "total current")]
    else:
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True)
        axes_to_plot = [
            (axes[0], total, "total current"),
            (axes[1], components[:, 0, :], "A current"),
            (axes[2], components[:, 1, :], "B current"),
        ]

    for ax, values, ylabel in axes_to_plot:
        for index in selected:
            ax.plot(time, values[:, index], linewidth=1.4, label=f"x={x_values[index]:.4g}")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    axes_to_plot[-1][0].set_xlabel("time")
    fig.suptitle(title)
    fig.tight_layout()
    output = Path(output_path)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved current time-trace plot to {output.resolve()}")


def plot_current_spacetime_map_from_npz(
    input_path: str | Path = "length_current_distribution.npz",
    output_path: str | Path = "current_spacetime_total.png",
    subtract_initial: bool = False,
) -> None:
    """Load a current npz result and plot the total-current spacetime map."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "current_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"current result is missing keys: {missing}")

    plot_current_spacetime_map(
        result["time_grid"],
        result["x_grid"],
        result["current_total"],
        output_path=output_path,
        title=f"Total current from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_current_time_slices_from_npz(
    input_path: str | Path = "length_current_distribution.npz",
    output_path: str | Path = "current_time_slices.png",
    time_indices: Iterable[int] | None = None,
    time_points: Iterable[float] | None = None,
    subtract_initial: bool = False,
) -> None:
    """Load a current npz result and plot total-current profiles at selected times."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "current_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"current result is missing keys: {missing}")

    plot_current_time_slices(
        result["time_grid"],
        result["x_grid"],
        result["current_total"],
        output_path=output_path,
        time_indices=time_indices,
        time_points=time_points,
        title=f"Current profiles from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_current_time_traces_from_npz(
    input_path: str | Path = "length_current_distribution.npz",
    output_path: str | Path = "current_time_traces.png",
    x_indices: Iterable[int] | None = None,
    x_points: Iterable[float] | None = None,
    subtract_initial: bool = False,
) -> None:
    """Load a current npz result and plot time traces at selected x points."""
    result = load_npz_result(input_path)
    required = ("time_grid", "x_grid", "current_total")
    missing = [key for key in required if key not in result]
    if missing:
        raise KeyError(f"current result is missing keys: {missing}")

    plot_current_time_traces(
        result["time_grid"],
        result["x_grid"],
        result["current_total"],
        current_components=result.get("current_components"),
        output_path=output_path,
        x_indices=x_indices,
        x_points=x_points,
        title=f"Current time traces from {Path(input_path).name}",
        subtract_initial=subtract_initial,
    )


def plot_current_standard_views_from_npz(
    input_path: str | Path = "length_current_distribution.npz",
    output_prefix: str | Path = "current",
    subtract_initial: bool = False,
) -> None:
    """Create the standard static current views: spacetime map, time slices, and time traces."""
    prefix = Path(output_prefix)
    plot_current_spacetime_map_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_spacetime_total.png"),
        subtract_initial=subtract_initial,
    )
    plot_current_time_slices_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_time_slices.png"),
        subtract_initial=subtract_initial,
    )
    plot_current_time_traces_from_npz(
        input_path,
        output_path=prefix.with_name(f"{prefix.name}_time_traces.png"),
        subtract_initial=subtract_initial,
    )
