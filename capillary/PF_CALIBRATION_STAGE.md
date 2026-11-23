# PF-SPH surface-tension, wetting, and bridge-resolution calibration

This is an **experimental** staged validation, not a claim that a capillary bridge
has been validated against Young–Laplace. Native simulations require your
SPHinXsys build (on macOS/CLion). The C++ solver continues to use the existing
Whebell-style PF forces; no new pressure or wetting model was substituted.

## Why a previously repulsive bridge may be possible

With a poorly wetting solid, the interfacial and pressure terms can give an
outward load. The previous `dx=0.05, t=5` bridge has a late-time left-body force
about `-0.05245` (negative = repulsive), with +0.2239 from direct PF interaction
and -0.2762 from other forces. The sign alone does not determine the equilibrium
contact angle, and we cannot infer a capillary contact angle from the two-sphere
force until `sff` and `sfs` are calibrated.

## Measurement methodology

1. **Free droplets:** measure conservative core-region virial pressure (including
   pressure/EOS and pair-force contributions). Use at least three *relaxed*
   droplet radii, fit `P_core = P_offset + gamma * (2/R_vol)`. Gas pressure is
   zero by this model's gauge. Do not fit only one radius; do not set the
   intercept to zero as a matter of convenience. The solver appends
   `core_eos_pressure`, `core_mean_density` and `core_particle_count` to its
   existing `pf_history.csv` for independent audits. The virial pressure
   estimator is unchanged. Require a populated interior core `R-2H` with
   `H=4*dx`. A paper-specific comparison is `gamma ~ 30.96*sff`, but the
   scripts never substitute this empirical relation for the measured slope.
2. **Sessile droplets:** for the same `dx`, `sff` and fluid parameters, sweep
   `sfs/sff`. For each completed run, reconstruct the fluid density isosurface,
   fit a volume-constrained zero-gravity spherical cap, and inspect late-time
   angle drift, cutoff sensitivity, density error, Mach number, and volume
   mismatch. The scripts will NOT automatically accept an unstable angle fit.
3. **Fixed-sphere bridge:** require a calibration JSON for the current `dx`
   before running. The bridge is initialized at a nominal fixed liquid radius
   and mass; verify the actual mass in the output when refining.
4. **Refinement:** repeat **both calibrations** at `dx=0.025`, using the coarse
   measured gamma as a *target*, not forcing equal raw interaction strengths.
   Measure `gamma` and angle again. Only compare bridge forces if gamma, angle,
   mass, and physical sphere geometry are consistent.

The numerical thresholds are intentionally conservative screening criteria, not
statistical confidence intervals or a substitute for a third grid level.

## Merge on your Mac and build

From your CLion project directory:

```bash
unzip -oq ~/Downloads/sphinxsys_capillary_calibration_full_source.zip -d /tmp/sph-calibration
rsync -av /tmp/sph-calibration/sphinxsys-capillary/ "$PWD/"
# Ensure Ninja recompiles new adapter, regardless of ZIP timestamps:
touch src/pf_validation.cpp
cmake --build cmake-build-release-vcpkgsph --target pf_capillary_validation -j 6
ctest --test-dir cmake-build-release-vcpkgsph -R '^pf_' --output-on-failure
python3 -m pytest tests/test_pf_calibration.py
```

Run commands from the project directory. All use the same executable:

```bash
EXE="$PWD/cmake-build-release-vcpkgsph/pf_capillary_validation"
```

## First: gamma at dx=0.05

```bash
python3 tools/pf_calibrate.py droplets \
  --executable "$EXE" --output results-pf-gamma-dx005 \
  --dx 0.05 --radii 0.65 0.8 0.95 --sff 0.001 \
  --end-time 10 --threads 4 --snapshot-interval 1
```

Inspect `results-pf-gamma-dx005/gamma_report.json`. If it reports
`not_converged`, inspect each droplet log and history, and extend the
relaxation (new output directory); do not proceed with the default coefficients.
A positive slope, small intercept and small residual are required.

## Second: target theta = 60 degrees at dx=0.05

```bash
python3 tools/pf_calibrate.py angles \
  --executable "$EXE" --output results-pf-angle-dx005 \
  --dx 0.05 --gamma-report results-pf-gamma-dx005/gamma_report.json \
  --sfs-ratios 1 1.5 2 3 4 --target-angle 60 \
  --end-time 12 --threads 4 --snapshot-interval 1
```

This generates a detailed `angle_report.json`. Only if a stable angle lies
within 5 degrees of the desired target does the script emit
`calibration.json`. If it emits `needs_refinement.json` instead, adjust the
`--sfs-ratios` list to bracket the target or extend equilibration.

## Third: calibrated fixed bridge at dx=0.05

```bash
python3 tools/pf_calibrate.py bridge \
  --executable "$EXE" --output results-pf-bridge-cal-dx005 \
  --dx 0.05 --calibration results-pf-angle-dx005/calibration.json \
  --sphere-radius 0.30 --gap 0.20 --liquid-radius 0.35 \
  --end-time 5 --threads 4 --snapshot-interval 0.5
```

## Fourth: redo gamma and angle at dx=0.025

```bash
COARSE_GAMMA=$(python3 -c 'import json;print(json.load(open("results-pf-gamma-dx005/gamma_report.json"))["gamma"])')

python3 tools/pf_calibrate.py droplets \
  --executable "$EXE" --output results-pf-gamma-dx0025 \
  --dx 0.025 --radii 0.40 0.50 0.60 --sff 0.001 \
  --target-gamma "$COARSE_GAMMA" --end-time 10 --threads 4
```

**Important:** If `gamma_report.json` says `target_mismatch`, rerun the finer
droplets at its `suggested_next_sff` (in a **new output directory**), and
use that successful gamma report for the subsequent sessile runs. This
suggestion assumes approximately linear gamma vs `sff` and must be verified.

Once finer gamma matches:

```bash
python3 tools/pf_calibrate.py angles \
  --executable "$EXE" --output results-pf-angle-dx0025 \
  --dx 0.025 --gamma-report results-pf-gamma-dx0025/gamma_report.json \
  --sfs-ratios 1 1.5 2 3 4 --target-angle 60 \
  --end-time 12 --threads 4

python3 tools/pf_calibrate.py bridge \
  --executable "$EXE" --output results-pf-bridge-cal-dx0025 \
  --dx 0.025 --calibration results-pf-angle-dx0025/calibration.json \
  --sphere-radius 0.30 --gap 0.20 --liquid-radius 0.35 \
  --end-time 5 --threads 4 --snapshot-interval 0.5
```

If finer gamma needs retuning, change the `--gamma-report` in the fine-angle
command to the **new** accepted report. Likewise, use the angle calibration
associated with that report.

## Compare two calibrated bridge resolutions

```bash
python3 tools/pf_calibrate.py compare \
  --coarse results-pf-bridge-cal-dx005/bridge \
  --fine results-pf-bridge-cal-dx0025/bridge \
  --coarse-cal results-pf-angle-dx005/calibration.json \
  --fine-cal results-pf-angle-dx0025/calibration.json \
  --output bridge_convergence.json
```

The comparison refuses mismatched surface tensions (>10%), contact angles
(>5 degrees), mass (>4%), and sphere geometry. Reported force agreement at
two resolutions is only a preliminary convergence indication, not a proven
asymptotic convergence order. The comparison also attempts the repository's
axisymmetric Young–Laplace reference for contact angles below 90 degrees,
reporting an unresolved reference rather than inventing one if its solution
cannot be found. A third resolution and a resolved independent equilibrium
benchmark remain necessary.

**Archives:** Each simulation stage writes a compressed `results-*.zip` next
to its output directory (without huge snapshots/checkpoints). Keep the original
folders for restart, and share the ZIP with the assistant for detailed review.

**Units:** Every calibration value is in the solver's current simulation units
(`rho0=1`, etc.). These scripts do not convert to SI, and `theta=60` is a
user-changeable target, not a claim about a measured material.
