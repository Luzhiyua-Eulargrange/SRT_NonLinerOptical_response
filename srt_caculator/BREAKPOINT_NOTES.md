# Breakpoint Notes

Date: 2026-06-26

## Current State

The project now has two separate RDM current entry points:

- `main_velocity_gauge.py`: velocity-gauge RDM current calculation.
- `main_length_gauge.py`: length-gauge RDM current calculation.

The default `main.py` calls the length-gauge entry point.

The shared lower-level modules are in place:

- `RDM_Velocity_Gauge.py`: velocity-gauge RDM evolution.
- `RDM_Length_Gauge.py`: length-gauge RDM evolution with Berry connection and covariant derivative.
- `Current.py`: current calculation from band-basis RDM trajectories.
- `Geometry.py`: velocity matrix, Berry connection, covariant derivative.
- `Debug_Tools.py`: band plot, current plot, RDM sanity checks.

The output design is:

- `*_band_result.npz`
- `*_band_structure.png`
- `*_rdm_current_result.npz`
- `*_current_convergence.npz`
- `*_current_evolution.png`

## Current Difficulty

Both velocity-gauge and length-gauge current convergence checks currently fail when comparing raw current traces between coarse and fine k grids.

Observed example:

```text
Length-gauge current convergence relative max difference: about 1.1
Velocity-gauge current convergence relative max difference: about 1.1
```

This suggests that the present convergence metric may not be physically appropriate. The current comparison is done directly on raw `J(t)`.

A likely issue is that raw current includes a k-grid-sensitive equilibrium or offset contribution. For optical response, it may be more meaningful to compare the response current:

```text
Delta J(t) = J(t) - J(0)
```

The next check should compare:

```python
relative_current_difference(
    fine_current - fine_current[0],
    coarse_current - coarse_current[0],
)
```

instead of comparing raw currents directly.

## Important Scientific Caution

The current operator for length gauge is currently evaluated using the same band velocity matrix path as the velocity-gauge current:

```text
J(t) = -e sum_k w_k Re Tr[v(k) rho(k,t)]
```

This is a reasonable first implementation, but it should be checked carefully against the paper's gauge-specific current expressions. The RDM evolution equations have been separated by gauge, but the current extraction still needs theoretical review.

## Next Suggested Steps

1. Replace the raw-current convergence check with a response-current convergence check using `J(t) - J(0)`.
2. Re-run `main_length_gauge.py` first.
3. If length gauge converges, save it as the primary reliable workflow.
4. Then run `main_velocity_gauge.py` and document whether it fails as expected due to velocity-gauge truncation sensitivity.
5. Review the length-gauge current operator against the paper before treating the numerical current as final physics.

## Files Most Likely To Touch Next

- `main_length_gauge.py`
- `main_velocity_gauge.py`
- `Current.py`
- `README.md`
