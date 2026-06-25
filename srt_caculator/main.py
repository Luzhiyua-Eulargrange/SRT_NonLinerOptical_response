import numpy as np

from Band_Solver import save_band_results, solve_bands
from Debug_Tools import debug_rdm_trajectory, plot_bands_and_print_eigenvectors, print_default_params
from RDM_Length_Gauge import propagate_length_gauge_rdm
from RDM_Velocity_Gauge import propagate_velocity_gauge_rdm
from config import make_k_grid, normalize_params


def main():
    print_default_params()

    params = normalize_params({
        "L": 2,
        "E0": 0.02,
        "omega": 0.5,
        "pulse_duration": 40.0,
        "num_time_steps": 201,
        "fermi_energy": 0.0,
    })

    k_grid, k_weight = make_k_grid(params, num_k=51)

    energies, eigenvectors = solve_bands(k_grid, params)
    print("Band calculation finished")
    print("k_grid shape:", k_grid.shape)
    print("energies shape:", energies.shape)
    print("eigenvectors shape:", eigenvectors.shape)

    plot_bands_and_print_eigenvectors(
        k_grid,
        energies,
        eigenvectors,
        output_path="band_structure.png",
        k_indices=(0,),
    )

    save_band_results(
        "band_result.npz",
        k_grid=k_grid,
        k_weight=k_weight,
        energies=energies,
        eigenvectors=eigenvectors,
    )
    print("Saved band results to band_result.npz")

    rdm_debug_params = dict(params)
    rdm_debug_params.update({
        "pulse_duration": 4.0,
        "t_switch": 1.0,
        "rtol": 1e-5,
        "atol": 1e-7,
    })
    rdm_debug_k_grid, _ = make_k_grid(rdm_debug_params, num_k=5)
    rdm_debug_energies, rdm_debug_eigenvectors = solve_bands(rdm_debug_k_grid, rdm_debug_params)
    debug_time_grid = np.linspace(0.0, rdm_debug_params["pulse_duration"], 5)

    time_grid_velocity, rho_velocity = propagate_velocity_gauge_rdm(
        rdm_debug_params,
        rdm_debug_k_grid,
        rdm_debug_energies,
        rdm_debug_eigenvectors,
        time_grid=debug_time_grid,
    )
    debug_rdm_trajectory(
        time_grid_velocity,
        rho_velocity,
        title="Velocity-gauge RDM debug check",
    )

    time_grid_length, rho_length = propagate_length_gauge_rdm(
        rdm_debug_params,
        rdm_debug_k_grid,
        rdm_debug_energies,
        rdm_debug_eigenvectors,
        time_grid=debug_time_grid,
    )
    debug_rdm_trajectory(
        time_grid_length,
        rho_length,
        title="Length-gauge RDM debug check",
    )


if __name__ == "__main__":
    main()
