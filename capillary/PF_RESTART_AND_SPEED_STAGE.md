# PF-SPH continuation: restart, reciprocal parallelism and performance audit

This source bundle is an **update to** the `sphinxsys-capillary` project,
not a new simulator or a completed moving-DEM coupling. It preserves the
Whebell et al. (2025) liquid-only PF-SPH equations and the pinned SPHinXsys
neighbor/storage infrastructure. The source contains no SPHinXsys checkout;
reuse your existing `deps/SPHinXsys` and local dependencies.

## What changed

- `src/pf_liquid.h`: only one pair distance calculation for pressure, projected
  viscosity and PF cohesion; invariant `H^-4` Wendland prefactor is evaluated
  once per interaction pass. The midpoint integrator no longer copies the large
  initial rate arrays when they are already available.
- `src/pf_parallel.h`: deterministic, mutually exclusive **TBB** pair shards.
  Each shard accumulates complete reciprocal pair contributions into local
  arrays, then shards are reduced in a fixed order. No pairwise atomics, data
  races, single-sided force approximations or momentum projections. The memory
  overhead is O(number_of_threads * number_of_samples). Native `--threads 1`
  retains serial accumulation. Speedup is not yet measured.
- `src/pf_restart.h`: versioned full-state checkpoints include the case and
  material/geometry settings, sample positions/velocities/densities/masses,
  global step and time, final-step dt and the **integrated linear and angular
  impulses on the wall**. The file is atomically replaced by rename (POSIX).
  The native executable checks particle counts, mass/type and wall geometry on
  restart. It recomputes neighbor lists and forces from the recovered state.
- `src/pf_snapshot.h`: backward-compatible recovery from existing
  `pf_snapshot_<step>.csv` and sibling `pf_snapshots.csv`; no need to start
  the 70-minute old run over. **Old CSVs do not contain wall impulse history.**
  The conservation ledger is rebased to zero at the resumed snapshot, and this
  fact is recorded in the run metadata. This is NOT an exact historical impulse
  reconstruction, though the CSV fluid state is restored.
- `src/pf_validation.cpp`: `--checkpoint-interval T`, `--resume FILE` and
  `--resume-snapshot CSV` options. New `pf_perf.csv` logs cumulative wall time,
  neighbor-graph time, force-evaluation time, call counts and mean milliseconds.
- `tools/run_pf_resume.py`: one short recovery run with fresh output. It never
  launches a two-case sweep and does not overwrite the archived trajectory.

## Rebuild (Linux)

From your existing `/media/egor/Data/sphinxsys-capillary` checkout, after
copying the archive's `sphinxsys-capillary/` **contents** over the matching
project files:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target pf_capillary_validation pf_restart_test pf_snapshot_test pf_liquid_test -j 6
ctest --test-dir build -R '^pf_' --output-on-failure
```

The native physics executable cannot be built by the standalone-only CMake mode.
If a native build fails, send the complete compiler error; do not start a long
physical validation until it builds and the native smoke tests pass.

## First run: short snapshot recovery only

The uploaded stopped trajectory reached step 26100 at t=7.81679, but its last
complete saved CSV snapshot is `pf_snapshot_23377.csv` at t=7.0002326908245482.
Recover from t=7.00023 and run to t=7.05 (about 0.05 additional simulation units):

```bash
python3 tools/run_pf_resume.py \
  --executable "$PWD/build/pf_capillary_validation" \
  --snapshot results-pf-wetting-dx004/pf_snapshot_23377.csv \
  --end-time 7.05 --threads 6 \
  --output results-pf-dx004-restart-smoke
```

The script writes `results-pf-dx004-restart-smoke.zip`, including the full
`pf_checkpoint.dat`, `pf_history.csv`, `pf_perf.csv`, and `run.log`.
Use a *different* `--output` for subsequent runs. The script's `--sound-speed 20
--dt-scale 1` defaults are **inferred** from the archived dt vs speed history
rather than independently recorded command metadata. The other defaults are
rho0=1, radius=0.6, dx=0.04, sff=0.001, sfs=0.0035, viscosity=0.05.
If your original dx=0.04 command used different settings, supply the actual
values explicitly before interpreting the continued trajectory.

To restart this **new** full-state checkpoint (with its saved wall impulse ledger),
start a separate new output directory and use:

```bash
mkdir -p results-pf-dx004-continued
cd results-pf-dx004-continued
../build/pf_capillary_validation \
  --case sessile --dx .04 --radius .6 --sff .001 --sfs .0035 \
  --viscosity .05 --sound-speed 20 --dt-scale 1 \
  --end-time 7.25 --threads 6 --record-every 100 --dump-every 0 \
  --checkpoint-interval .1 \
  --resume ../results-pf-dx004-restart-smoke/pf_checkpoint.dat
```

All resumed output is **a new segment**: histories and snapshot indices do not
silently concatenate into the old run. Full checkpoints preserve integrated wall
impulse from the start of their own ledger; for the legacy-recovered segment,
that ledger begins at the t=7 snapshot.

## Validation status / gating

Run locally in this environment: four standalone C++ tests, including serial
vs sharded accumulation and bit-for-bit fixed-graph trajectory split/restart,
and 21 Python regression tests. A compile-only TBB API shim also checked the
`pf_parallel.h` template and sharded numerical equivalence; it is NOT a real
oneTBB or SPHinXsys build.

**NOT RUN HERE:** full native SPHinXsys build, native CTest smoke, wall-clock
thread-scaling benchmark, full dx=0.04 trajectory, or moving DEM validation.
Do not claim a measured speedup or quantitative wetting convergence from the
standalone tests.

Next checks after the short native smoke: compare `--threads 1` and `--threads 6`
for the same recovered t=7.00 to t=7.05 segment, inspect `pf_perf.csv` to see
whether force evaluation or the native neighbor rebuild dominates, and verify
mass, density, angular/linear residuals, apparent contact-angle drift. Then
proceed to moving curved-solid surface samples with equal/opposite pressure,
viscous and PF loads, DEM mass/inertia, DEM contacts and capillary bridge
formation/rupture tests. This bundle does **not** implement those stages.
