# Local verification, 2026-11-02

- Linux x86-64, GCC 13.3.0, C++17, Release.
- Clean pinned SPHinXsys v1.0-beta.08 library built, then capillary_validation linked.
- C++ syntax check passed against actual dependency headers.
- Three independent reference tests passed: exact cylindrical bridge, geometric
  similarity of Young–Laplace bridge, spherical-cap shape/angle reconstruction.
- CTest execution check passed.
- Four cases (Laplace, sessile, fixed bridge, moving bridge) each completed TWO
  acoustic steps at coarse dx=0.3. Runner analysis/archive path completed.
- Liquid mass unchanged in these short tests. Volume and momentum are separately
  recorded; no conservation projection or mass correction was applied.

No equilibrium benchmark was run locally. These coarse tests deliberately do
not meet the resolution/equilibration requirements for force or angle accuracy.
Nonzero supported momentum residuals exist in the short outputs; they are not
suppressed, and this project does not claim a validated capillary reaction or
production DEM coupling. Run the spatial/time/sound-speed study and review the
actual results before using particle loads. Raw short evidence is under verification/.


## Focused contact-angle update (2026-11-03)

Current update: eight Python tests pass (analytic caps at 60/90/120 degrees with
interior layers, degenerate input, angle history, snapshot-index integration,
two-case runner selection, and prior Young–Laplace reference tests).
The actual `time_window.h` helper compiles with GCC/C++17 and passes the standalone
endpoint/non-progress regression for end times .01, .3, 20 and 40. These are clock
operations only, not fluid trajectories. Four archived two-step outputs also
postprocess successfully with the updated analyzer.

The full modified SPHinXsys adapter has NOT been rebuilt or run in this environment:
CMake, Eigen, TBB and Simbody build dependencies are absent. The earlier build logs
refer to the original adapter, not this update. Rebuild locally and run the focused
`--smoke` check before long trajectories. No new force or wetting physics is claimed.


## Force-accounting diagnostic update (2026-11-03)

Ten Python tests and the standalone C++ clock regression pass. The Python suite
checks signed external-CSF accounting without erasing the raw momentum defect and
supports old result files. Dry-run command generation selects exactly one case for
the requested force audit. No fluid trajectory was run here. The full newly
instrumented adapter has not been rebuilt in this dependency-free environment.
The uploaded t=40 runs verify the preceding adapter version, not this instrumentation.


## Acoustic viscosity update (2026-11-03)

- Ten Python reference/analysis/runner tests pass.
- Both standalone C++17 regressions compile and pass: endpoint handling, and
  viscous row coefficients plus unequal-mass pair momentum conservation,
  energy dissipation, bounded velocities and timestep convergence.
- The full updated `src/validation.cpp` passes GCC C++17 `-fsyntax-only` against
  the pinned upstream headers, Eigen 3.4, TBB 2021.11, Simbody 3.7 and Boost 1.83
  headers. Headers were extracted to a temporary prefix; no system installation
  was needed. This checks template/API compatibility but is not a linked build.
- Runner command/archive plumbing passes using a fake adapter and archived
  smoke data. It does not constitute a new SPH run.
- No fluid trajectory or long benchmark was run locally. The requested t=2
  acoustic-viscosity audit must still establish whether gas spikes are reduced.
- The wetting-reaction/momentum-conservation issue remains unresolved.


## Capillary energy derivation/reference update (2026-11-03)

Fifteen Python tests pass, including five new checks of an independent geometric
energy-gradient oracle: coordinate finite differences, Euclidean invariance,
force/torque identities, rigid virtual work, spherical-cap uplift and fixed-volume
Young equilibrium. No production capillary operator changed. This verifies the
reference force algebra, not a new SPH surface-area discretization or moving DEM
coupling. No fluid trajectory was run or requested for this stage.


## Independent gas-free PF-SPH pair module (2026-11-03)

The actual header compiles with GCC, C++17, -O2 -DNDEBUG -Wall -Wextra -pedantic.
Its standalone test passes profile, potential/force gradient, configuration energy
gradient, Euclidean covariance, force/torque reciprocity, rigid virtual work,
unequal-mass impulses, scaling and invalid-input checks. Checks remain enabled
under NDEBUG. Existing fifteen Python tests also pass. No fluid trajectory was run.
The production adapter was not changed to gas-free; the CMake pair-only option is
provided for an isolated dependency-free build. Direct compiler execution was tested;
CMake configure/build for this option was not run in this environment.

## Liquid-only adapter stage (2026-11-03)

- Added `pf_liquid.h` (published WCSPH/central viscosity/continuity/midpoint core)
  and `pf_validation.cpp` (gas-free native SPHinXsys adapter).
- Recompiled both new C++ tests with g++ C++17, O2, NDEBUG, Wall/Wextra/pedantic;
  both ran successfully. Liquid checks cover pressure–continuity power balance,
  nonpositive viscous power, reciprocal forces/torques, integrated wall impulses,
  fixed wall/mass invariants and a tiny two-step 3D liquid lattice.
- C++17 syntax/template check of `pf_validation.cpp` against pinned SPHinXsys and
  extracted Eigen/TBB/Simbody headers passed. See `verification/pf_adapter_syntax.log`.
- Existing Python regressions: 15 passed in 1.884 s.
- New two-case Python runner: byte-compilation and dry-run command generation passed.
- Full linking and native SPHinXsys smoke execution were NOT performed here:
  installed CMake and full dependency libraries are unavailable. Two native CTests
  are provided for the user's build. No long trajectories were run locally.
- No physical calibration or moving-DEM validation is claimed. See PF_LIQUID_RUNS.md.

## Density-isosurface diagnostic stage

- 19 Python tests passed in 1.547 s (4 new independent interface checks).
- Processed the existing uploaded wetting final frame at grid/dx=.25 and .125;
  these are offline measurements, not new simulations. Reports are included in
  verification/pf_wetting_interface_grid_*.json.
- Current C++ snapshot I/O block compiled with C++17/Wall/Wextra/pedantic and ran
  successfully with a two-frame CSV index. Full SPHinXsys rebuild was not performed
  for the added output flag/snapshot changes; old syntax logs describe the prior
  adapter revision. No long computational benchmarks were run.

## Focused wetting runner stage

- Added two-run relaxation/compressibility runner; c10 scale=.5 and c20 scale=1
  use matched acoustic bounds and identical geometry/strengths at t=15.
- Added observed angle drift summary with geometry screening and at least three
  late-time valid frames; equilibrium remains unestablished.
- 21 Python regressions passed in 1.565 s. Dry-run command generation, Python
  byte-compilation and git diff --check passed. No C++ changes this stage and no
  new fluid simulations were run locally.
