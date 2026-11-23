# Liquid-only PF-SPH pilot: build and run

The new executable is **`pf_capillary_validation`**. The older
`capillary_validation` executable remains the archived two-phase model.
There are no gas particles in the new executable. Do not use the old runner for it.

From the project root, with your existing SPHinXsys dependencies installed:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target pf_capillary_validation pf_liquid_test pf_capillary_test -j 6
ctest --test-dir build -R '^pf_' --output-on-failure
```

The two native adapter smoke tests run at most two steps each. They also compare
native neighbor coverage against an exhaustive pair search on their tiny initial
systems. They are executable checks, not physical validation.

Then run exactly two longer pilot trajectories:

```bash
python3 tools/run_pf_validation.py \
  --executable "$PWD/build/pf_capillary_validation" \
  --output results-pf-first --end-time 5 --threads 6
```

Send back **`results-pf-first.zip`**. The runner leaves VTP files locally in each
case's `output/` directory and excludes those bulky files from the diagnostic zip.
The archive contains histories, final liquid coordinates, run settings and logs.
Use a new output name on each attempt. Add `--dump-every 0` to disable VTP.
`--dry-run` prints commands without creating files or launching simulations.
NumPy/SciPy are optional for the final convex-envelope sphere/angle fit; the
simulation and basic CSV diagnostics do not depend on Python scientific packages.

The cases are a free spherical droplet (R=.6, dx=.04) and a sessile initial
hemisphere (R=.6, dx=.06), in the dimensionless units of the earlier adapter.
The sessile substrate is flat and stationary. Both have zero gravity and zero
initial velocity. Material parameters are rho0=1, artificial sound speed=10,
dynamic viscosity=.05, s_ff=.001 and s_fs=.001. These last two values are pilot
strengths, NOT a claimed surface tension/contact angle calibration. Vary `--sff`
and `--sfs` explicitly; there is intentionally no `--theta` input yet.
The native command supports `--radius`, `--dx`, `--dt-scale`, `--record-every`,
`--dump-every` and `--steps`; see `pf_capillary_validation --help`.

## What is implemented

Whebell et al., JFM 1019 A37 (2025), DOI 10.1017/jfm.2025.10601,
equations 2.4–2.21, with the already-tested pair potential:

- Fixed liquid sample masses; density evolves by the continuity equation. No
  surface density summation reset, global mass/volume restoration, gas, CSF or
  prescribed interface-normal rotation is used.
- Tait EOS with exponent 7. Negative EOS pressure is retained, as in the paper.
- Quintic Wendland kernel W=21/(2*pi*H^3)*(1-q)^4*(1+4q), H=4dx.
  Native SPHinXsys neighbor support is explicitly enlarged to this radius before
  constructing relations. Forces always use current coordinates, not cached
  native distances or quadrature weights.
- Unique pair pressure force
  F_ij^P= -m_i*m_j/(rho_i*rho_j)*(P_i+P_j)*grad_i W.
- Central projected viscosity
  F_ij^mu=10*mu*m_i*m_j/(rho_i*rho_j)*(v_ij dot e_ij)/r*W'(r)*e_ij.
  W'<=0 makes pair viscous power nonpositive. Unlike a componentwise velocity
  drag, this operator conserves pair angular momentum.
- Fluid-solid dummy samples enter the SAME pressure/viscosity/PF summations.
  Wall samples have fixed density rho0 and velocity zero. Every force has its
  opposite accumulated on that wall sample, including capillarity.
- The published modified position-Verlet/continuity integration (2.15–2.16).
  Neighbors refresh at old and midpoint configurations. The final current rates
  are reused at the start of the next step, avoiding redundant force evaluation.
  Midpoint wall impulses are accumulated over the SAME completed timestep.
- Paper acoustic/acceleration/viscous timestep restrictions, strengthened by an
  advective speed term and explicit projected-viscosity row-sum limit. No
  velocity/force/density clipping or momentum projection is applied.

The source splits the independently tested core (`src/pf_liquid.h`) from the
SPHinXsys adapter (`src/pf_validation.cpp`). SPHinXsys supplies storage, lattice
sampling, cell lists, neighbor relations, VTP I/O and TBB infrastructure.
**The new reciprocal pair loop currently runs serially.** `--threads` controls
native TBB operations; it does not make the whole solver scale across that many
cores. Parallelizing reciprocal accumulation is deferred until this baseline
passes fluid tests. There is no OpenMP dependency or upstream library patch.

## What to inspect

`pf_history.csv` records fixed mass, compressible volume sum(m/rho), kinetic and
pair energies, density errors, speed, reciprocal wall loads and system momentum
errors. Mass conservation does NOT imply volume is constant in WCSPH. Pressure
work is consistent with the continuity operator at the differential level;
finite-timestep total energy conservation is not asserted. Pair energy alone is
not total energy: compressibility has an internal-energy contribution as well.

The sessile momentum check is fluid momentum plus integrated force on the fixed
wall. The angular check includes wall reaction moments about the global origin.
A fixed substrate is an external support, not a freely moving DEM body.

For the free droplet, the reported virial pressure includes BOTH pressure and PF
forces, excludes viscosity, and follows the paper's directed-pair half-sum at
r_v=R_volume-2H. It is reported as NaN when that core is unresolved. The Python
summary reports a one-radius apparent sigma=mean(P_virial*R_volume/2), for
calibration only. The paper's sigma=30.96*s_ff mapping is printed as a comparison,
not used as an input or enforced. A multi-radius fit and spatial virial plateau
are still required for Laplace validation. Coarse pilot results may have substantial
interface/support-size bias. No mechanical surface-tension claim follows from
merely completing this run.

The optional final angle fit uses the existing convex-envelope estimator. It
assumes one approximately spherical convex droplet; it is NOT suitable for
separated bridges or arbitrary aggregates. Wetting has a finite interaction
range and a numerical wall layer, so sensitivity to interface location, wall
position and resolution must be assessed before assigning a calibrated angle.
The runner never declares equilibrium automatically.

## Checks and scope

Executed locally: dependency-free C++ pair and liquid-core regression tests,
including a two-step 3D liquid-only lattice, pressure–continuity power balance,
viscous dissipation, wall impulse accounting, linear and angular momentum.
The new full adapter passes C++17 syntax/template checking against the pinned
SPHinXsys headers. All existing 15 Python regressions pass.

**The full SPHinXsys executable was not linked or run in this environment** because
the complete installed build dependencies/CMake are unavailable. The short native
CTest runs above are the first check on your machine. No long trajectory was run
here. This is an implemented, unvalidated liquid-only solver baseline, not a
validated DEM coupling. Moving solids, overlaps, curved-solid wetting calibration,
DEM contacts and joining/breaking bridge benchmarks are subsequent stages.

## Density-isosurface angle measurement (new)

Install the postprocessing dependencies and analyze an existing result without
rerunning its fluid simulation or rebuilding:

```bash
python3 -m pip install -r requirements.txt
python3 tools/pf_interface.py results-pf-wetting --grid-dx 0.125 --mesh
```

This writes `pf_interface_analysis.json` and optionally the fluid-only triangle
mesh `pf_interface.vtp`, which ParaView can open directly. It reconstructs
sum_j m_j W(x-x_j)=rho0/2 with the same H=4dx Wendland kernel. Only liquid samples
are included. It does not change fluid mass, density, forces or trajectories.
Use `--grid-dx 0.1` for the paper's grid spacing; the default 0.25 is a faster pilot.
A memory guard rejects grids over 30 million nodes. Reconstruction is offline.

The zero-gravity Young–Laplace solution is a spherical cap. The fit uses the known
volume sum(m/rho) to link curvature and angle, fits an independent vertical
interface shift, and uses inverse cylindrical-radius weighting. It does not
choose the shift to obtain a target angle. Data within 1, 1.5 and 2 sample spacings
of the physical wall are excluded in separate fits to expose near-wall bias.
The 1.5dx fit is reported as primary. This exclusion is an additional diagnostic
choice, not a claim to exactly reproduce every detail of the paper's fitter.
Finite-gravity profiles, nonconvex drops and bridges are not supported.

The report includes the reconstructed enclosed volume, its difference from
sum(m/rho), connected-surface count, fit residuals and cutoff sensitivity. More
than one surface component or a volume mismatch above 10% makes it ineligible
for calibration. Passing these geometry screens never establishes equilibrium
or proves physical validity. An underresolved wall layer can remain biased.

Your saved stronger-adhesion run (sff=.001, sfs=.0035, t=5) gives:

| Grid spacing / dx | Fitted angle | Isosurface volume error |
| --- | --- | --- |
| 0.25 | 61.451 degrees | -6.145% |
| 0.125 | 61.478 degrees | -6.096% |

The reconstruction-grid effect is small. This measurement does NOT recover the
paper's roughly 45-degree calibration, so measurement choice alone does not
explain that discrepancy. The simulation was still relaxing; fluid spatial
resolution, compressibility, initialization and boundary discretization remain
unresolved sensitivities. Do not retune adhesion to hide these errors.

For FUTURE native runs, rebuild `pf_capillary_validation` and add
`--snapshot-interval 1` (physical simulation time). This records `pf_snapshot_*.csv`
and `pf_snapshots.csv` independently of VTP cadence; zero disables it. Frame times
are actual completed timesteps and are not forced onto an exact output schedule.
The two-case runner enables this interval by default. Afterward `pf_interface.py`
also writes `pf_angle_history.csv` from those snapshots. A historical run with only
a final CSV cannot provide a measured angle history retroactively.

The native adapter now explicitly marks liquid state changes for VTP output.
Without this notification the library could emit only the initial frame despite
later dump requests. This output fix does not affect numerical dynamics.

Verification: 19 Python checks pass, including kernel normalization, an analytic
single-kernel isosurface volume, and exact shifted spherical caps at 40, 75 and
120 degrees. The added snapshot I/O block was compiled and run separately in
C++17. The complete adapter was not rebuilt here this turn because its installed
header/library dependencies are unavailable; earlier full-adapter syntax checks
predate this output-only change. No new fluid trajectories were run locally.

## Focused wetting relaxation/compressibility checks

The next runner launches exactly TWO sessile runs, both with dx=.06, initial
R=.6, sff=.001, sfs=.0035, mu=.05 and end time 15:

| Case | Sound speed | Timestep scale | Purpose |
| --- | --- | --- | --- |
| c10_dt_half | 10 | 0.5 | Longer relaxation; compare its t≈5 frame against the earlier c=10, scale=1 result |
| c20_dt_full | 20 | 1.0 | Stiffer EOS at approximately the same acoustic timestep |

The acoustic bound .15*dx/c times timestep scale is equal for the two cases.
The velocity, acceleration or viscous restrictions can still make actual steps
slightly different. This is a controlled sensitivity check, not proof of temporal
or incompressible convergence. All forces, particle masses and viscosity remain
unchanged. Grid refinement and bridge validation follow after these results.

After extracting the updated sources into the project, rebuild the snapshot-enabled
executable and run the short native checks:

```bash
python3 -m pip install -r requirements.txt
cmake --build build --target pf_capillary_validation pf_capillary_test pf_liquid_test -j 6
ctest --test-dir build -R '^pf_' --output-on-failure
```

Then launch the two trajectories:

```bash
python3 tools/run_pf_wetting_checks.py \
  --executable "$PWD/build/pf_capillary_validation" \
  --output results-pf-wetting-checks --threads 6
```

The script checks Python postprocessing dependencies before launching, records
command lines and an executable hash, runs each trajectory, reconstructs the
saved interfaces and writes angle histories and late-time drift summaries. It
never declares equilibrium automatically. Invalid geometry frames are excluded
from the angle-drift calculation but remain in the history. At least three valid
frames in the final quarter are required to report a slope.

Default snapshot interval is 1 simulation time unit and reconstruction grid/dx=.25.
There is no VTP output by default; optional `--dump-every 2000` enables native
particle output. The diagnostic archive contains CSV/JSON/logs, including the
liquid snapshots needed for independent reconstruction. Existing simulation output
is preserved: select a new output name for a fresh pair of runs.

Send back `results-pf-wetting-checks.zip`. No free-droplet rerun or parameter sweep
is launched. To see the two commands without running anything, add `--dry-run`.
To repeat offline analysis at finer reconstruction spacing without rerunning fluid:

```bash
python3 tools/run_pf_wetting_checks.py \
  --output results-pf-wetting-checks --analyze-only --grid-dx 0.125
```

Analysis-only refreshes the reports and diagnostic archive but preserves raw
simulation files. A failed simulation or fit retains partial output and is included
in the archive with its error; the runner exits nonzero after packaging it.

Local verification for this runner update: 21 Python tests passed in 1.565 s,
including the paired parameter contract and late-time drift on a known synthetic
angle history with an invalid geometry frame. Dry-run and byte-compilation passed.
No fluid trajectories, long benchmarks or C++ physics changes were made this stage.
