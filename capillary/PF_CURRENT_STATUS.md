# Calibration review — 2026-11-05

Update 2026-11-06: general N-particle DEM development now continues with provisional
parameters; see [PF_DEM_COUPLING.md](PF_DEM_COUPLING.md). The calibration findings
below remain unchanged and do not prevent running development examples.

Continue from the neighbor-acceleration full source, which supersedes the earlier
restart/parallelization package. This review changes Python diagnostics only;
it does not modify native forces, integration, neighbor lists, or checkpoints.

## Verified from the other chat's uploaded results

- Moving-sphere example completed to t=0.2. Linear and angular residuals were
  about 3.3e-18 and 8.7e-19. This is a numerical coupling check, not physical validation.
- Neighbor benchmark: native wall time 3.218 s, Verlet wall time 0.522 s (6.17x).
  Solver-only speedup was 6.69x; graph builds dropped from 422 to 2.
  Coordinates, density and velocity agreed within approximately 4.4e-15.
  These timings apply to that short benchmark, not every aggregate.
- Checkpoints, fixed/free two-sphere bridges, rigid-body quaternions, and cached
  neighbors are already present in this baseline. Do not redo those stages.
- The dx=0.05, sff=0.001 three-radius droplet runs completed to t=10, but did
  not pass the surface-tension calibration. The fixed bridge test therefore
  remains uncalibrated; it is not evidence of a validated capillary-force law.

## What the pressure review found

The existing fit gives gamma=0.0466798 and pressure intercept=-0.0251232.
The intercept is 26.3% of the mean absolute measured pressure, above the existing
20% gate. The paper's 30.96*sff comparison is 0.03096; it is not a replacement
for an actual calibration at this discretization.

The largest droplet's late pressure standard deviation is 29.6% of its mean.
A small linear slope alone had classified this oscillating signal as stationary.
The updated diagnostic additionally requires standard deviation / absolute mean
<=5%, and reports explicit failed checks. This is a conservative numerical screen,
not a statistical confidence interval or proof of equilibrium. Correlated samples
and slow oscillations still require longer-window comparisons.

An independent reconstruction from each final particle CSV reproduces the native
stored final virial pressure within 1.5e-15. EOS and PF contributions have opposite
signs, each around magnitude 4, while the net pressure is about 0.06–0.12.
The resulting cancellation factors are about 69, 86 and 133. This makes small
sampling/discretization changes significant relative to the net pressure.
Agreement verifies this calculation, not the physical validity of its estimator.

The native estimator currently uses the single core R_volume-2H. The paper's
radial averaging for sufficiently large droplets should be evaluated separately.
The offline tool scans core radii and reports geometric versus sampled volume
sensitivity. It DOES NOT substitute a volume renormalization or tune the force law
to force agreement with the paper. Only one final frame is available per uploaded
case; the scan cannot establish a time-averaged gamma or resolve temporal oscillations.

## Reproduce the review without running any simulations

From the repository root, with the gamma result directory extracted:

```bash
python3 -m pytest -q tests/test_pf_calibration.py tests/test_pf_pressure_audit.py
python3 tools/pf_calibrate.py droplets --dx .05 --sff .001 \
  --radii .65 .8 .95 --output results-pf-gamma-dx005 --analyze-only
python3 tools/pf_pressure_audit.py results-pf-gamma-dx005 \
  --output results-pf-pressure-audit --plot
```

The test command needs pytest; the offline plot needs matplotlib. Numerical
postprocessing dependencies remain listed in requirements.txt. Neither analysis
command launches an executable. Analyze-only refreshes JSON analyses; an existing
ZIP is preserved, so collect the refreshed JSON separately if sending new results.

The included `verification/other_chat_review/` contains the reanalysis, final-frame
pressure profiles, plot, and the 12-test log from this review. No native trajectory
or long benchmark was run locally. Native compilation was not attempted here.

## Next decision

Resolve pressure sampling and temporal averaging before selecting sfs from contact
angle or validating bridge force. First inspect the existing largest-droplet history
and sample pressure at multiple core radii through time; do not repeat an entire
angle sweep. Additional runs, if needed, should target that unresolved measurement
and retain checkpoints and frame histories. The compact uploaded gamma archive
contains final frames but not the full checkpoints needed for exact continuation;
those may remain in the original local result directories.
