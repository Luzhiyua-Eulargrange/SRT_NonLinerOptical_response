"""Velocity-gauge RDM current calculation.

This entry point is intentionally explicit. Velocity-gauge current convergence
is known to be fragile for truncated band spaces, so failed convergence raises
an error instead of silently saving a questionable result.
"""

from __future__ import annotations

import numpy as np

from Band_Solver import save_band_results, solve_bands
from Current import save_rdm_current_results, total_current_from_band_rdm
from Debug_Tools import debug_rdm_trajectory, plot_bands_and_print_eigenvectors, plot_current_evolution, print_default_params
from RDM_Velocity_Gauge import propagate_velocity_gauge_rdm
from config import make_k_grid, normalize_params


def relative_current_difference(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return max-norm relative difference between two current traces."""
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    if reference.shape != candidate.shape:
        raise ValueError("current traces must have the same shape")
    scale = max(float(np.max(np.abs(reference))), 1e-12)
    return float(np.max(np.abs(candidate - reference)) / scale)


def solve_velocity_current(params: dict, num_k: int, time_grid: np.ndarray) -> dict:
    """Solve velocity-gauge RDM and current on one k grid."""
    run_params = dict(params)
    run_params["num_k_rdm"] = int(num_k)
    k_grid, k_weight = make_k_grid(run_params)
    energies, eigenvectors = solve_bands(k_grid, run_params)
    time_grid_out, rho_trajectory = propagate_velocity_gauge_rdm(
        run_params,
        k_grid,
        energies,
        eigenvectors,
        time_grid=time_grid,
    )
    current = total_current_from_band_rdm(
        run_params,
        k_grid,
        k_weight,
        rho_trajectory,
        eigenvectors,
    )
    return {
        "k_grid": k_grid,
        "k_weight": k_weight,
        "energies": energies,
        "eigenvectors": eigenvectors,
        "time_grid": time_grid_out,
        "rho_trajectory": rho_trajectory,
        "current": current,
    }


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
    convergence_tol = 0.20

    band_k_grid, band_k_weight = make_k_grid(params, num_k=params["num_k_band"])
    band_energies, band_eigenvectors = solve_bands(band_k_grid, params)
    plot_bands_and_print_eigenvectors(
        band_k_grid,
        band_energies,
        band_eigenvectors,
        output_path="velocity_band_structure.png",
        print_eigenvectors=False,
    )
    save_band_results(
        "velocity_band_result.npz",
        k_grid=band_k_grid,
        k_weight=band_k_weight,
        energies=band_energies,
        eigenvectors=band_eigenvectors,
    )

    time_grid = np.linspace(0.0, params["pulse_duration"], params["num_time_steps"])
    coarse_num_k = max(5, params["num_k_rdm"] // 2)
    fine_num_k = params["num_k_rdm"]

    print(f"Solving velocity-gauge coarse grid Nk={coarse_num_k}")
    coarse = solve_velocity_current(params, coarse_num_k, time_grid)
    debug_rdm_trajectory(coarse["time_grid"], coarse["rho_trajectory"], "Velocity coarse RDM")

    print(f"Solving velocity-gauge fine grid Nk={fine_num_k}")
    fine = solve_velocity_current(params, fine_num_k, time_grid)
    debug_rdm_trajectory(fine["time_grid"], fine["rho_trajectory"], "Velocity fine RDM")

    convergence_error = relative_current_difference(fine["current"], coarse["current"])
    print("Velocity-gauge current convergence relative max difference:", convergence_error)
    if convergence_error > convergence_tol:
        raise RuntimeError(
            "Velocity-gauge current did not converge: "
            f"relative max difference {convergence_error:.6g} exceeds {convergence_tol:.6g}. "
            "This is expected for some truncated-band velocity-gauge calculations; "
            "increase num_k_rdm or compare with length gauge."
        )

    save_rdm_current_results(
        "velocity_rdm_current_result.npz",
        "velocity",
        fine["k_grid"],
        fine["k_weight"],
        fine["time_grid"],
        fine["rho_trajectory"],
        fine["current"],
        fine["energies"],
    )
    np.savez(
        "velocity_current_convergence.npz",
        coarse_num_k=coarse_num_k,
        fine_num_k=fine_num_k,
        time_grid=fine["time_grid"],
        coarse_current=coarse["current"],
        fine_current=fine["current"],
        relative_max_difference=convergence_error,
        tolerance=convergence_tol,
    )
    plot_current_evolution(fine["time_grid"], fine["current"], "velocity_current_evolution.png")


if __name__ == "__main__":
    main()
