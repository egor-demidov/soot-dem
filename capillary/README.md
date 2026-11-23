For installation into the original soot-dem project, use the additive installer:
[SOOT_BUILD_MERGE.md](SOOT_BUILD_MERGE.md). Do not overwrite soot-dem's root
CMakeLists.txt with this component's standalone build file.

CMake builds the standalone PF-SPH/DEM targets by default. No SPHinXsys bootstrap
is required. To build the older SPHinXsys validation adapters explicitly, enable
`-DPF_BUILD_SPHINXSYS_VALIDATION=ON` and set `SPHINXSYS_SOURCE_DIR` to an existing
checkout. The legacy `CAPILLARY_PAIR_TESTS_ONLY=ON` flag remains supported.



Actual soot-dem neck/contact coupling: [SOOT_DEM_COUPLING.md](SOOT_DEM_COUPLING.md).

General moving-particle coupling and ready-to-run XML examples: [PF_DEM_COUPLING.md](PF_DEM_COUPLING.md).

Current baseline and recovered calibration results: [PF_CURRENT_STATUS.md](PF_CURRENT_STATUS.md).
### Surface tension and contact angle calibration (experimental)

See [PF_CALIBRATION_STAGE.md](PF_CALIBRATION_STAGE.md) for the staged
multiradius Laplace-pressure, sessile contact-angle, and two-resolution
fixed-bridge verification workflow. The native solver now appends core EOS
pressure/density/population fields to `pf_history.csv` for pressure audits.
The Python workflow uses `tools/pf_calibrate.py` and does not claim physical
validation until native droplet and wetting measurements have been obtained.

## Latest experimental two-sphere stage (November 2026)

`pf_capillary_validation --case bridge` adds a finite liquid bridge between
two fixed or moving DEM spheres, force decomposition and version-3 restart.
See [PF_TWO_SPHERE_BRIDGE_STAGE.md](PF_TWO_SPHERE_BRIDGE_STAGE.md) for the
macOS/CLion commands, test runner and important physical limitations.

# Experimental moving-sphere DEM–SPH update

See [PF_MOVING_DEM_STAGE.md](PF_MOVING_DEM_STAGE.md). This stage introduces
an independently tested reciprocal 6DOF rigid sphere core, a native (not yet
runtime-verified) `--case immersed` adapter, and v2 full-state DEM checkpointing.
**This is not yet a validated two-sphere capillary-bridge implementation.**

# Latest update: PF-SPH checkpoint, threading and legacy snapshot recovery

See [PF_RESTART_AND_SPEED_STAGE.md](PF_RESTART_AND_SPEED_STAGE.md) for the new
`--checkpoint-interval`, `--resume`, `--resume-snapshot`, `pf_perf.csv`, and
one-run dx=0.04 recovery workflow. Native SPHinXsys performance and moving-DEM
physics are not yet validated. The earlier baseline instructions follow below.

# New: gas-free PF-SPH pilot

The new liquid-only executable is `pf_capillary_validation`. Build/run instructions
and scope are in [PF_LIQUID_RUNS.md](PF_LIQUID_RUNS.md). The old
`capillary_validation` executable and its runners still use the two-phase model.

# SPHinXsys capillary validation — independent clean CMake project

This project tests whether an existing SPH implementation is suitable for resolved
capillary bridges around moving DEM spheres. It does not depend on your soot-dem
repository, waLBerla, or our previous mixture kernels. No external DEM integration
is claimed. SPHinXsys is pinned, and its physics source is unmodified.

## Build

Dependencies: C++17 compiler, CMake >=3.20, Git, Eigen3, oneTBB, Simbody and Boost
(headers and program_options). Python diagnostics require NumPy and SciPy.

Ubuntu/Debian:
```bash
sudo apt-get install cmake ninja-build g++ libeigen3-dev libtbb-dev \
  libsimbody-dev libboost-dev libboost-program-options-dev
```
macOS/Homebrew:
```bash
brew install cmake ninja eigen tbb simbody boost
```
The local build check is Linux/GCC. Apple Silicon has not been tested here;
the old dependency may need portability work, which must be tracked explicitly.

From this project directory:
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python tools/bootstrap.py
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target capillary_validation -j 4
python -m unittest discover -s tests -v
python tools/run_validation.py --smoke --threads 2
```
On macOS, if Homebrew packages are not found, add
`-DCMAKE_PREFIX_PATH="$(brew --prefix)"` to the configure command. No OpenMP
installation is needed: the selected upstream path uses TBB.

## Physical validation runs (run on your machine)

For the next stage, run **two contact-angle cases only**, at dx=0.06 and
60/120 degrees. Rebuild after applying this update:

```bash
cmake --build build --target capillary_validation -j 4
python tools/run_contact_angle.py --executable build/capillary_validation --threads 6 --end-time 40
```

This starts fresh trajectories (restart is not implemented). It writes no VTP
files, records scalar diagnostics every 100 acoustic steps, and saves liquid-only
CSV positions about once per simulation time unit plus initial/final states.
Postprocessing generates `angle_history.csv`, final-quarter angle slope/range,
and separate liquid/gas Mach and density diagnostics. Send `results-ca-long.zip`.
To run just one case, add `--angles 120`; to extend duration use `--end-time 80`
and a fresh `--output` directory. The time is not an equilibrium guarantee.

A cheap adapter/output check before the long run is:

```bash
python tools/run_contact_angle.py --executable build/capillary_validation --smoke --output results-ca-smoke
```

The older `run_validation.py --refinement` sweep remains available for later
convergence studies; do not rerun it at this stage.

1. **Laplace pressure:** spherical droplet R=0.6 in gas, zero gravity. Reference
   delta-p=2 sigma/R. The final independent surface fit estimates R.
2. **Wetting:** an initially hemispherical droplet on z=0 relaxes toward input
   contact angle 60/90/120 degrees. A geometric sphere fit to convex-envelope vertices above 2dx estimates the angle
   independently of imposed normals. It assumes a single approximately convex drop;
   cutoff and radial-offset sensitivity are reported, not calibrated away.
3. **Fixed bridge:** two spheres R=0.6, surface separation 0.3, initial liquid
   cylinder radius 0.4, offset x=0.04. Independent axisymmetric Young–Laplace
   shooting solves a symmetric equilibrium branch for measured liquid volume.
   Compare the upper-sphere axial force against that reference only after the
   offset bridge has relaxed back toward axisymmetry. The reference can fail to
   find a branch and must then be reported as unresolved, never a validation pass.
4. **Moving bridge:** the same pair is released. Uses upstream pressure and
   viscous loads, finite dry mass and spherical inertia; records force and torque,
   motion, and supported momentum. This is an experimental coupling audit. It stops with an explicit failure if spheres reach contact: a DEM contact law is not yet integrated.

All quantities use a consistent nondimensional system: rho_liquid=1, rho_gas=.001,
mu_liquid=.05, mu_gas=.00005, sigma=.008. Artificial sound speed=10. Do not map
these to your nanometre aggregate until the dimensionless tests pass. Wetting
case and pressure case require explicit gas particles in this implementation.

## Viscous time integration update

The default `--viscosity-update acoustic` recomputes both phase viscous accelerations
before every pressure substep. The original schedule is retained explicitly as
`--viscosity-update advection` for reproducibility; it does not apply the new limit.
Surface-tension geometry/forcing still updates at the advection frequency. Force
accumulators reset before recomputation, then the stored CSF force is restored;
viscosity is never added repeatedly to an unreset accumulator.

For the pinned scalar viscosity operator, each velocity component obeys
`du_i/dt = sum_j b_ij (u_j-u_i)`, with nonnegative coefficients. The numerator
`rho_i*b_ij` is `-2*mu*dW_ij*V_j/(r_ij+0.01*h)` for inner links and
`-4*mu_effective*dW_ij*V_j/(r_ij+0.01*h)` for phase/mirrored-wall links.
The phase viscosity is the upstream harmonic mean. The acoustic timestep also
satisfies `dt <= 0.5/max_i(sum_j b_ij)`. This gives a convex-combination bound for
an isolated explicit viscous update. It is not a proof of stability for the
coupled pressure, interface, transport-correction and moving-boundary algorithms.
Cached numerators use exactly the neighbor data held fixed by the SPH operators;
current densities enter the limit at each acoustic step. No velocity clipping,
additional viscosity or global momentum projection is introduced.

Compare against your existing t=2 audit with one new short run:

```bash
cmake --build build --target capillary_validation -j 4
python tools/run_contact_angle.py --executable build/capillary_validation \
  --angles 60 --end-time 2 --audit-every 100 --viscosity-update acoustic \
  --threads 6 --output results-viscosity-audit
```

Send `results-viscosity-audit.zip`. The run metadata records the chosen schedule;
scalar history includes the measured rate, viscous limit and last acoustic step.
The previous wetting-reaction gap remains unresolved intentionally. First determine
whether the gas velocity outliers disappear with consistent viscous timestepping.

## Next diagnostic: one short force audit

After rebuilding, run only one 60-degree trajectory to t=2:

```bash
python tools/run_contact_angle.py --executable build/capillary_validation \
  --angles 60 --end-time 2 --audit-every 100 --output results-force-audit --threads 6
```

Send `results-force-audit.zip`. This is a diagnostic run, not an equilibrium test.
The audit instrumentation records:

- CSF linear and angular impulses integrated each acoustic step (all runs).
- With `--audit-every N`: momenta before/after each liquid/gas pressure and density
  update, total prior force, and the pre-update pressure/viscous wall load.
- At scalar-record times: the eight fastest gas particles, positions, velocities,
  densities, accelerations, solid signed distance, nearest liquid-center distance,
  and neighbor counts cached at the most recent relation rebuild.

`analysis.json` includes both the original supported momentum defect and a separate
balance treating the recorded continuum surface force (CSF) as external. The latter
is **only an accounting test**; it does not establish closed-system conservation or
supply a physical wetting reaction to DEM. No force projection, velocity clipping,
artificial damping, or reaction-force correction has been added.

The pinned upstream wetting example also omits gas pressure feedback in the liquid
pressure operator while evolving gas with a multiphase operator. This is a deliberate
upstream choice, not a new adapter omission; it is another reason not to assume
exact phase-pair momentum conservation. Separate wall-load sampling precedes pressure
initialization and does not exactly match the gas extended wall operator. The audit
must isolate these contributions before selecting a conservative moving-solid model.

## Evidence and acceptance

CSV records separate liquid mass (sum m), liquid volume (sum m/rho), compressibility,
and momentum. Constant particle mass does not imply constant volume. The primary shape fit uses convex-envelope particle positions; the former thick
surface-indicator band fit is retained as a legacy diagnostic. Radial offsets
0, half-dx and dx are sensitivity tests, not calibrated corrections or rigorous
uncertainty bounds. Inspect the shape: a convex envelope cannot resolve concavities
or establish that there are no detached fragments.
The offset must become irrelevant under spatial refinement. Small bridges and
near-contact surfaces can remain unresolved even when the solver runs.

Suggested initial accuracy targets, to review rather than enforce blindly:
- relative liquid mass drift <1e-12 and volume variation <1%;
- max density deviation <2%, Mach <0.1;
- equilibrated Laplace pressure and bridge force error <5%;
- equilibrated contact angle error <5 degrees;
- errors decrease under spatial refinement and are insensitive to increased
  artificial sound speed. For moving coupling, additionally halve time steps
  before interpreting a momentum defect as a missing force.

Examine the last 20% of history for equilibrium: pressure, angle/shape, load and
kinetic energy must stop drifting. No automatic acceptance is based merely on
elapsed time or a zero exit code. `analysis.json` always says NOT_VALIDATED until
these comparisons are reviewed. See UPSTREAM_AUDIT.md before using particle loads.

The outer container is fixed. Its load impulse is included as support_impulse;
compare P_fluid+P_dry+support_impulse with its initial value, and similarly for
angular momentum. In a fixed bridge, both spheres are also external supports.
Only the moving case includes dry-body momenta. Pressure and viscous loads are
sampled before the fluid substep, so this is a diagnostic balance, not an exact
algebraic identity. A small CSF sum alone does not establish particle-force accuracy.

Outputs: history.csv, loads.csv, run.json, analysis.json, liquid snapshots in CSV
and standard SPHinXsys VTP particle output readable in ParaView. VTP represents
particles, not a reconstructed fluid-only triangle mesh. The runner includes a
source/executable hash manifest. No long benchmark is run automatically by CMake.

## Git

Keep this as a separate repository while evaluating:
```bash
git init
git add CMakeLists.txt README.md UPSTREAM_AUDIT.md src tools tests requirements.txt .gitignore
git commit -m "Add pinned SPHinXsys capillary validation project"
```
The dependency is ignored and reproducibly fetched by bootstrap.py; its exact
commit is tracked in bootstrap and CMake. A submodule can replace this later.

SPHinXsys upstream: https://github.com/Xiangyu-Hu/SPHinXsys (Apache-2.0).


## Capillary coupling derivation and reference checks

See `CAPILLARY_COUPLING_DESIGN.md` for the energy-based force/torque framework and
its remaining SPH discretization requirements. `research/capillary_energy_reference.py`
is a geometric oracle for force, torque and virtual work; it is not imported by the
simulation and does not add a triangle interface tracker. Run its small tests with
`python -m unittest discover -s tests -v`. No new fluid trajectory is requested for
this derivation update. Production capillarity remains unchanged.


## Gas-free PF-SPH module (first implementation stage)

`src/pf_capillary.h` implements the published Whebell et al. (2025) central pair
force profiles, their analytic potential, reciprocal accumulation and rigid-load
reduction. See `PF_SPH_IMPLEMENTATION.md` for units, assumptions and tests. Build
`pf_capillary_test`, or configure with `-DCAPILLARY_PAIR_TESTS_ONLY=ON` to avoid all
SPHinXsys dependencies. This module is not yet wired into `capillary_validation`;
that executable remains the existing two-phase adapter. Do not launch another long
fluid run for this update. The next step is a compatible liquid-only adapter and
surface-tension/wettability calibration, not another custom area estimator.

## Optional neighbor-list acceleration

See [PF_NEIGHBOR_ACCELERATION.md](PF_NEIGHBOR_ACCELERATION.md) for the
cell-linked and displacement-triggered Verlet backends, native fallback,
validation, and the macOS/CLion A/B benchmark runner. `verlet` is the default;
use `--neighbor-mode native` to reproduce the original neighbor search.

## Original soot-dem XML inputs

Updated original droplet and levelset inputs are in `inputs/`. See
[XML_INPUT_COMPATIBILITY.md](XML_INPUT_COMPATIBILITY.md) for units, includes,
commands, output controls, and the explicit provisional SPH calibration settings.

For the fixed-boundary penetration investigation and two fresh controls, see
[FIXED_BOUNDARY_AUDIT.md](FIXED_BOUNDARY_AUDIT.md).

Optional initial fluid-only relaxation (particles frozen): see [FLUID_EQUILIBRATION.md](FLUID_EQUILIBRATION.md).
