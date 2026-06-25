import numpy as np

from config import normalize_params, make_k_grid
from Band_Solver import solve_bands
from Current import total_current
from Debug_Tools import print_default_params


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

    time_grid, current = total_current(
          params,
          k_grid=k_grid,
          k_weight=k_weight,
      )

    print("Current calculation finished")
    print("time_grid shape:", time_grid.shape)
    print("current shape:", current.shape)

    np.savez(
        "result.npz",
        k_grid=k_grid,
        energies=energies,
        time_grid=time_grid,
        current=current,
    )

    print("Saved result to result.npz")


if __name__ == "__main__":
    main()
