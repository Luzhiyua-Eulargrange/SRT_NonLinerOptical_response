"""Length-gauge RDM convergence calculation."""

from __future__ import annotations

import numpy as np

from Band_Solver import save_band_results, solve_bands
from Current_Distribution import (
    current_convergence as real_current_convergence,
    current_result_from_rdm_result,
    save_current_result,
)
from Density import (
    density_convergence as real_density_convergence,
    density_result_from_rdm_result,
    save_density_result,
)
from Debug_Tools import debug_rdm_trajectory, plot_bands_and_print_eigenvectors, print_default_params
from Debug_Tools import plot_current_spacetime_map, plot_current_time_slices, plot_current_time_traces
from Debug_Tools import plot_density_spacetime_map, plot_density_time_slices, plot_density_time_traces
from RDM_Length_Gauge import propagate_length_gauge_rdm
from config import make_k_grid, normalize_params


def relative_max_difference(reference: np.ndarray, candidate: np.ndarray, floor: float = 1e-12) -> float:
    """Return max-norm relative difference between two arrays."""
    reference = np.asarray(reference)
    candidate = np.asarray(candidate)
    if reference.shape != candidate.shape:
        raise ValueError("arrays must have the same shape")
    scale = max(float(np.max(np.abs(reference))), float(floor))
    return float(np.max(np.abs(candidate - reference)) / scale)


def density_matrix_purity(rho_trajectory: np.ndarray) -> np.ndarray:
    """Return Tr(rho^2) for every time and k point."""
    rho = np.asarray(rho_trajectory, dtype=np.complex128)
    if rho.ndim != 4 or rho.shape[-1] != rho.shape[-2]:
        raise ValueError("rho_trajectory must have shape (Nt, Nk, Nb, Nb)")
    return np.einsum("tkij,tkji->tk", rho, rho).real


def density_matrix_eigenvalues(rho_trajectory: np.ndarray) -> np.ndarray:
    """Return sorted Hermitian eigenvalues of rho for every time and k point."""
    rho = np.asarray(rho_trajectory, dtype=np.complex128)
    if rho.ndim != 4 or rho.shape[-1] != rho.shape[-2]:
        raise ValueError("rho_trajectory must have shape (Nt, Nk, Nb, Nb)")

    nt, nk, nb, _ = rho.shape
    values = np.empty((nt, nk, nb), dtype=float)
    for it in range(nt):
        for ik in range(nk):
            values[it, ik] = np.linalg.eigvalsh(rho[it, ik])
    return values


def nested_fine_indices(coarse_k_grid: np.ndarray, fine_k_grid: np.ndarray) -> np.ndarray:
    """Return fine-grid indices matching a nested endpoint-excluding coarse grid."""
    coarse_k_grid = np.asarray(coarse_k_grid, dtype=float)
    fine_k_grid = np.asarray(fine_k_grid, dtype=float)
    if coarse_k_grid.ndim != 1 or fine_k_grid.ndim != 1:
        raise ValueError("k grids must be one-dimensional")
    if fine_k_grid.size % coarse_k_grid.size != 0:
        raise ValueError("fine k grid size must be an integer multiple of coarse k grid size")

    stride = fine_k_grid.size // coarse_k_grid.size
    indices = np.arange(coarse_k_grid.size) * stride
    if not np.allclose(fine_k_grid[indices], coarse_k_grid, rtol=0.0, atol=1e-12):
        raise ValueError("fine k grid is not nested with the coarse k grid")
    return indices


def density_matrix_convergence(coarse: dict, fine: dict) -> dict[str, float]:
    """
    Compare RDM trajectories on nested k grids.

    The full rho matrix comparison assumes both runs use the same band gauge
    convention. Diagonal populations, traces, purity, and rho eigenvalues are
    included as more robust diagnostics.
    """
    if not np.allclose(coarse["time_grid"], fine["time_grid"], rtol=0.0, atol=1e-12):
        raise ValueError("coarse and fine runs must use the same time grid")

    fine_indices = nested_fine_indices(coarse["k_grid"], fine["k_grid"])
    coarse_rho = np.asarray(coarse["rho_trajectory"], dtype=np.complex128)
    fine_rho = np.asarray(fine["rho_trajectory"], dtype=np.complex128)[:, fine_indices]
    if coarse_rho.shape != fine_rho.shape:
        raise ValueError("coarse rho and restricted fine rho have inconsistent shapes")

    coarse_diag = np.diagonal(coarse_rho, axis1=-2, axis2=-1).real
    fine_diag = np.diagonal(fine_rho, axis1=-2, axis2=-1).real
    coarse_trace = np.trace(coarse_rho, axis1=-2, axis2=-1)
    fine_trace = np.trace(fine_rho, axis1=-2, axis2=-1)
    coarse_purity = density_matrix_purity(coarse_rho)
    fine_purity = density_matrix_purity(fine_rho)
    coarse_eigenvalues = density_matrix_eigenvalues(coarse_rho)
    fine_eigenvalues = density_matrix_eigenvalues(fine_rho)

    return {
        "rho_relative_max_difference": relative_max_difference(fine_rho, coarse_rho),
        "rho_max_abs_difference": float(np.max(np.abs(coarse_rho - fine_rho))),
        "diagonal_relative_max_difference": relative_max_difference(fine_diag, coarse_diag),
        "diagonal_max_abs_difference": float(np.max(np.abs(coarse_diag - fine_diag))),
        "trace_max_abs_difference": float(np.max(np.abs(coarse_trace - fine_trace))),
        "purity_relative_max_difference": relative_max_difference(fine_purity, coarse_purity),
        "purity_max_abs_difference": float(np.max(np.abs(coarse_purity - fine_purity))),
        "rho_eigenvalue_relative_max_difference": relative_max_difference(fine_eigenvalues, coarse_eigenvalues),
        "rho_eigenvalue_max_abs_difference": float(np.max(np.abs(coarse_eigenvalues - fine_eigenvalues))),
    }


def print_density_matrix_convergence(metrics: dict[str, float]) -> None:
    """Print density-matrix convergence metrics in a stable order."""
    print("Length-gauge RDM convergence diagnostics:")
    for key in (
        "rho_relative_max_difference",
        "rho_max_abs_difference",
        "diagonal_relative_max_difference",
        "diagonal_max_abs_difference",
        "trace_max_abs_difference",
        "purity_relative_max_difference",
        "purity_max_abs_difference",
        "rho_eigenvalue_relative_max_difference",
        "rho_eigenvalue_max_abs_difference",
    ):
        print(f"  {key}: {metrics[key]}")


def print_real_density_convergence(metrics: dict[str, float]) -> None:
    """Print real-space density convergence metrics in a stable order."""
    print("Length-gauge real-space density convergence diagnostics:")
    for key in (
        "density_total_relative_max_difference",
        "density_total_max_abs_difference",
        "density_component_relative_max_difference",
        "density_component_max_abs_difference",
        "charge_per_cell_relative_max_difference",
        "charge_per_cell_max_abs_difference",
    ):
        print(f"  {key}: {metrics[key]}")


def print_real_current_convergence(metrics: dict[str, float]) -> None:
    """Print real-space current convergence metrics in a stable order."""
    print("Length-gauge real-space current convergence diagnostics:")
    for key in (
        "current_total_relative_max_difference",
        "current_total_max_abs_difference",
        "current_component_relative_max_difference",
        "current_component_max_abs_difference",
        "paramagnetic_component_relative_max_difference",
        "paramagnetic_component_max_abs_difference",
        "diamagnetic_component_relative_max_difference",
        "diamagnetic_component_max_abs_difference",
        "current_per_cell_relative_max_difference",
        "current_per_cell_max_abs_difference",
        "average_current_relative_max_difference",
        "average_current_max_abs_difference",
    ):
        print(f"  {key}: {metrics[key]}")


def solve_length_rdm(params: dict, num_k: int, time_grid: np.ndarray) -> dict:
    """Solve length-gauge RDM on one k grid."""
    run_params = dict(params)
    run_params["num_k_rdm"] = int(num_k)
    k_grid, k_weight = make_k_grid(run_params)
    energies, eigenvectors = solve_bands(k_grid, run_params)
    time_grid_out, rho_trajectory = propagate_length_gauge_rdm(
        run_params,
        k_grid,
        energies,
        eigenvectors,
        time_grid=time_grid,
    )
    return {
        "k_grid": k_grid,
        "k_weight": k_weight,
        "energies": energies,
        "eigenvectors": eigenvectors,
        "time_grid": time_grid_out,
        "rho_trajectory": rho_trajectory,
    }


def save_length_rdm_result(filename: str, result: dict) -> None:
    """Save a length-gauge RDM result."""
    np.savez(
        filename,
        gauge=np.asarray("length"),
        k_grid=np.asarray(result["k_grid"], dtype=float),
        k_weight=float(result["k_weight"]),
        time_grid=np.asarray(result["time_grid"], dtype=float),
        rho_trajectory=np.asarray(result["rho_trajectory"], dtype=np.complex128),
        energies=np.asarray(result["energies"], dtype=float),
        eigenvectors=np.asarray(result["eigenvectors"], dtype=np.complex128),
    )


def save_rdm_convergence(
    filename: str,
    coarse: dict,
    fine: dict,
    metrics: dict[str, float],
    tolerance: float,
) -> None:
    """Save nested-grid RDM convergence diagnostics."""
    np.savez(
        filename,
        coarse_num_k=coarse["k_grid"].size,
        fine_num_k=fine["k_grid"].size,
        coarse_k_grid=np.asarray(coarse["k_grid"], dtype=float),
        fine_k_grid=np.asarray(fine["k_grid"], dtype=float),
        time_grid=np.asarray(fine["time_grid"], dtype=float),
        coarse_rho_trajectory=np.asarray(coarse["rho_trajectory"], dtype=np.complex128),
        fine_rho_trajectory=np.asarray(fine["rho_trajectory"], dtype=np.complex128),
        tolerance=float(tolerance),
        **metrics,
    )


def save_real_density_convergence(
    filename: str,
    coarse_density: dict,
    fine_density: dict,
    metrics: dict[str, float],
    tolerance: float,
) -> None:
    """Save real-space density convergence diagnostics."""
    np.savez(
        filename,
        x_grid=np.asarray(fine_density["x_grid"], dtype=float),
        time_grid=np.asarray(fine_density["time_grid"], dtype=float),
        coarse_density_total=np.asarray(coarse_density["density_total"], dtype=float),
        fine_density_total=np.asarray(fine_density["density_total"], dtype=float),
        coarse_density_components=np.asarray(coarse_density["density_components"], dtype=float),
        fine_density_components=np.asarray(fine_density["density_components"], dtype=float),
        coarse_particles_per_cell=np.asarray(coarse_density["particles_per_cell"], dtype=float),
        fine_particles_per_cell=np.asarray(fine_density["particles_per_cell"], dtype=float),
        tolerance=float(tolerance),
        **metrics,
    )


def save_real_current_convergence(
    filename: str,
    coarse_current: dict,
    fine_current: dict,
    metrics: dict[str, float],
    tolerance: float,
) -> None:
    """Save real-space current convergence diagnostics."""
    np.savez(
        filename,
        x_grid=np.asarray(fine_current["x_grid"], dtype=float),
        time_grid=np.asarray(fine_current["time_grid"], dtype=float),
        vector_potential=np.asarray(fine_current["vector_potential"], dtype=float),
        coarse_current_total=np.asarray(coarse_current["current_total"], dtype=float),
        fine_current_total=np.asarray(fine_current["current_total"], dtype=float),
        coarse_current_components=np.asarray(coarse_current["current_components"], dtype=float),
        fine_current_components=np.asarray(fine_current["current_components"], dtype=float),
        coarse_paramagnetic_current_components=np.asarray(
            coarse_current["paramagnetic_current_components"], dtype=float
        ),
        fine_paramagnetic_current_components=np.asarray(
            fine_current["paramagnetic_current_components"], dtype=float
        ),
        coarse_diamagnetic_current_components=np.asarray(
            coarse_current["diamagnetic_current_components"], dtype=float
        ),
        fine_diamagnetic_current_components=np.asarray(
            fine_current["diamagnetic_current_components"], dtype=float
        ),
        coarse_current_per_cell=np.asarray(coarse_current["current_per_cell"], dtype=float),
        fine_current_per_cell=np.asarray(fine_current["current_per_cell"], dtype=float),
        coarse_average_current=np.asarray(coarse_current["average_current"], dtype=float),
        fine_average_current=np.asarray(fine_current["average_current"], dtype=float),
        tolerance=float(tolerance),
        **metrics,
    )


def main() -> None:
    print_default_params()
    params = normalize_params({
        "L": 2,
        "E0": 0.02,
        "omega": 0.5,
        "pulse_duration": 20.0,
        "t_switch": 4.0,
        "num_time_steps": 81,
        "num_k_band": 51,
        "num_k_rdm": 21,
        "fermi_energy": 0.0,
        "rtol": 1e-6,
        "atol": 1e-8,
    })
    convergence_tol = 0.01
    density_num_x = 201

    band_k_grid, band_k_weight = make_k_grid(params, num_k=params["num_k_band"])
    band_energies, band_eigenvectors = solve_bands(band_k_grid, params)
    plot_bands_and_print_eigenvectors(
        band_k_grid,
        band_energies,
        band_eigenvectors,
        output_path="length_band_structure.png",
        print_eigenvectors=False,
    )
    save_band_results(
        "length_band_result.npz",
        k_grid=band_k_grid,
        k_weight=band_k_weight,
        energies=band_energies,
        eigenvectors=band_eigenvectors,
    )

    time_grid = np.linspace(0.0, params["pulse_duration"], params["num_time_steps"])
    coarse_num_k = params["num_k_rdm"]
    fine_num_k = 2 * coarse_num_k

    print(f"Solving length-gauge coarse grid Nk={coarse_num_k}")
    coarse = solve_length_rdm(params, coarse_num_k, time_grid)
    debug_rdm_trajectory(coarse["time_grid"], coarse["rho_trajectory"], "Length coarse RDM")

    print(f"Solving length-gauge fine grid Nk={fine_num_k}")
    fine = solve_length_rdm(params, fine_num_k, time_grid)
    debug_rdm_trajectory(fine["time_grid"], fine["rho_trajectory"], "Length fine RDM")

    rdm_convergence_metrics = density_matrix_convergence(coarse, fine)
    print_density_matrix_convergence(rdm_convergence_metrics)

    print(f"Reconstructing unit-cell density on Nx={density_num_x}")
    coarse_density = density_result_from_rdm_result(params, coarse, num_x=density_num_x)
    fine_density = density_result_from_rdm_result(params, fine, num_x=density_num_x)
    density_convergence_metrics = real_density_convergence(fine_density, coarse_density)
    print_real_density_convergence(density_convergence_metrics)

    convergence_error = density_convergence_metrics["density_total_relative_max_difference"]
    if convergence_error > convergence_tol:
        print(
            "Length-gauge real-space density did not converge: "
            f"relative max difference {convergence_error:.6g} exceeds {convergence_tol:.6g}."
        )

    print(f"Reconstructing unit-cell current on Nx={density_num_x}")
    coarse_current = current_result_from_rdm_result(params, coarse, num_x=density_num_x)
    fine_current = current_result_from_rdm_result(params, fine, num_x=density_num_x)
    current_convergence_metrics = real_current_convergence(fine_current, coarse_current)
    print_real_current_convergence(current_convergence_metrics)

    current_convergence_error = current_convergence_metrics["current_total_relative_max_difference"]
    if current_convergence_error > convergence_tol:
        print(
            "Length-gauge real-space current did not converge: "
            f"relative max difference {current_convergence_error:.6g} exceeds {convergence_tol:.6g}."
        )

    save_length_rdm_result("length_rdm_result.npz", fine)
    save_density_result("length_density_result.npz", fine_density)
    save_current_result("length_current_distribution.npz", fine_current)
    plot_density_spacetime_map(
        fine_density["time_grid"],
        fine_density["x_grid"],
        fine_density["density_total"],
        output_path="density_spacetime_total.png",
        title="Length-gauge total density n(x,t)",
    )
    plot_density_spacetime_map(
        fine_density["time_grid"],
        fine_density["x_grid"],
        fine_density["density_total"],
        output_path="density_response_spacetime_total.png",
        title="Length-gauge density response delta n(x,t)",
        subtract_initial=True,
    )
    plot_density_time_slices(
        fine_density["time_grid"],
        fine_density["x_grid"],
        fine_density["density_total"],
        output_path="density_response_time_slices.png",
        title="Length-gauge density response profiles",
        subtract_initial=True,
    )
    plot_density_time_traces(
        fine_density["time_grid"],
        fine_density["x_grid"],
        fine_density["density_total"],
        density_components=fine_density["density_components"],
        output_path="density_response_time_traces.png",
        title="Length-gauge density response time traces",
        subtract_initial=True,
    )
    plot_current_spacetime_map(
        fine_current["time_grid"],
        fine_current["x_grid"],
        fine_current["current_total"],
        output_path="current_spacetime_total.png",
        title="Length-gauge total current j(x,t)",
    )
    plot_current_spacetime_map(
        fine_current["time_grid"],
        fine_current["x_grid"],
        fine_current["current_total"],
        output_path="current_response_spacetime_total.png",
        title="Length-gauge current response delta j(x,t)",
        subtract_initial=True,
    )
    plot_current_time_slices(
        fine_current["time_grid"],
        fine_current["x_grid"],
        fine_current["current_total"],
        output_path="current_response_time_slices.png",
        title="Length-gauge current response profiles",
        subtract_initial=True,
    )
    plot_current_time_traces(
        fine_current["time_grid"],
        fine_current["x_grid"],
        fine_current["current_total"],
        current_components=fine_current["current_components"],
        output_path="current_response_time_traces.png",
        title="Length-gauge current response time traces",
        subtract_initial=True,
    )
    save_rdm_convergence(
        "length_rho_convergence.npz",
        coarse,
        fine,
        rdm_convergence_metrics,
        convergence_tol,
    )
    save_real_density_convergence(
        "length_density_convergence.npz",
        coarse_density,
        fine_density,
        density_convergence_metrics,
        convergence_tol,
    )
    save_real_current_convergence(
        "length_current_convergence.npz",
        coarse_current,
        fine_current,
        current_convergence_metrics,
        convergence_tol,
    )


if __name__ == "__main__":
    main()
