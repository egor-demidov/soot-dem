# Density-failure restart investigation

This update adds read-only reporting; no force law, density integration, timestep
limit, checkpoint format, input signature, or stopping threshold changes.

Reinstall this component with the additive installer and rebuild **the executable
used by your run**, e.g. from your soot-dem root:

```bash
cmake --build cmake-build-release-gcc-12 --target pf_soot_dem -j 6
bash capillary/tools/run_density_restart.sh \
  capillary_levelset_small.xml \
  results-sph-levelset/dem_checkpoint.dat \
  "$PWD/cmake-build-release-gcc-12/pf_soot_dem" \
  results-density-debug
```

The script restarts for at most 5000 accepted steps, stopping at 3.4e-10 s.
It is intended for the checkpoint near 0.322 ns discussed in this run. For other
checkpoint times use run_pf_dem.py with the appropriate --end-time/--steps.
Choose a fresh output directory. Source checkpoint is not overwritten.
Exit 1 is expected if the existing operating screen is reproduced; diagnostics
are printed afterwards even when the solver exits unsuccessfully.

`density_extrema.csv` has two rows per ordinary logging event: minimum and maximum
liquid density. Columns include sample ID, rho/rho0, Tait EOS pressure, density
rate, position, velocity, nearest body index, surface gap and neighbor counts.
The gap is min_j(distance-to-center_j - radius_j), i.e. distance outside the sphere
union or a negative penetration indicator. Counts exclude self and use the
actual compact support; buried boundary sites are excluded. These are geometric
counts, not normalized kernel-support measures.

`continuity_divergence = -(drho/dt)/rho` uses the existing SPH continuity rate.
Positive means density depletion/expansion; negative means compression. It is
not an independent continuum velocity-gradient estimate. Time is seconds,
pressure Pa, gaps and coordinates m, divergence s^-1 with SI inputs.

`last_dt` is the previous accepted timestep, zero on the initial restart record.
It does not identify the active stability constraint. Sample IDs persist in
this formulation and allow tracking an extreme as it moves.

When the density/Mach operating screen trips, `failure_liquid.vtp` and
`failure_particles.vtp` contain that failing, completed integration state even
if it lies between regular dump events. They are visualization snapshots, not
checkpoints. Existing checkpoint remains the latest one saved at the configured
cadence. Color failure_liquid.vtp by Density to locate the extreme; compare its
position and the reported surface gap with particle geometry.

Send the script's printed final rows and preferably density_extrema.csv/run.log.
No need to upload the full trajectory. The new diagnostics establish the sign
and location of the density error; they do not establish a physical fix.
