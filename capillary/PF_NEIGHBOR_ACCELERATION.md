# PF-SPH neighbor-list acceleration (experimental native adapter)

## What changed

`pf_capillary_validation` now has three runtime neighbor backends:

- `--neighbor-mode native` retains the existing SPHinXsys contact/inner list rebuilding on each evaluation. Use this as a reference.
- `--neighbor-mode cells` uses a separate deterministic 3D hashed cell-linked search at the true cutoff (`H = 4 dx`), rebuilding each evaluation.
- `--neighbor-mode verlet` (new default) builds a candidate list with search cutoff `H + skin`. It reuses this list until **any** particle, including rigid DEM surface samples, has moved by at least `skin / 2` relative to the last list build.

The `--neighbor-skin-factor F` option defines `skin = F * H` (default 0.25, allowed 0 to 1). With `F = 0` the Verlet implementation rebuilds each time. All candidates, including those beyond H, are passed to the unchanged pair operator, which evaluates current distances and returns exactly zero when `r >= H`. No force laws, liquid mass, EOS, or integrator settings were changed. Pair lists are always sorted by particle index, without duplicates, for deterministic reduction order.

SPHinXsys particle fields are synchronized for VTP output, or during the legacy native neighbor search. There is no reason to rebuild its cell lists for a cached pair evaluation.

## Install in an existing CLion project

Extract the **patch** archive to a temporary directory and merge it with your existing project; do not replace the `deps/SPHinXsys` directory, `.venv`, or CMake build directory. From the CLion terminal at `/Users/egor/ClionProjects/dem-sph`:

```bash
unzip -oq ~/Downloads/sphinxsys_capillary_neighbor_acceleration_patch.zip -d /tmp/sph-neighbors
rsync -av /tmp/sph-neighbors/sphinxsys-capillary/ "$PWD/"
touch src/pf_validation.cpp CMakeLists.txt
cmake --build cmake-build-release-vcpkgsph --target pf_capillary_validation pf_neighbors_test pf_bridge_test -j 6
ctest --test-dir cmake-build-release-vcpkgsph -R 'pf_(neighbors_regression|bridge_regression|bridge_adapter_smoke)' --output-on-failure
```

`touch` is useful on macOS if an archive preserved timestamps older than the existing build. In CLion use **Tools → CMake → Reload CMake Project** if `pf_neighbors_test` does not appear immediately.

## Run a controlled A/B benchmark

```bash
python3 tools/benchmark_pf_neighbors.py \
  --executable "$PWD/cmake-build-release-vcpkgsph/pf_capillary_validation" \
  --output results-pf-neighbor-benchmark \
  --dx 0.05 --end-time 0.15 --threads 4 --skin-factor 0.25
```

This runs two **fresh** fixed-sphere bridges with identical physical parameters, one with `native` neighbors and one with `verlet`. It writes `results-pf-neighbor-benchmark.zip`, with separate histories, performance data, run logs, final liquid coordinates, and `neighbor_comparison.json`.

Benchmark requirements:

1. Both runs must reach the requested end time and the same integration step count.
2. Max absolute position, velocity and density discrepancies should be near the floating-point rounding scale (the script conservatively checks `1e-7` in the simulation units, and `1e-7` in left-sphere force).
3. The cached run should report many `neighbor_reuses` and relatively few `neighbor_builds`.
4. Compare **whole solver time** (`wall_seconds`) as well as neighbor search time; a larger skin gives fewer rebuilds but more candidate pairs, so maximum cache reuse does not necessarily imply fastest runtime.

Only after the native A/B check should you use the cached mode for calibration runs. To force old behavior, pass `--neighbor-mode native` to the executable. The calibration scripts use the new default unless explicitly modified.

## Tuning

Try skin factors `0.10`, `0.20`, `0.25`, and `0.40` with a **different** `--output` directory each time. Prefer the factor with the lowest overall solver wall time, subject to matching final states. Factors are fractions of the support radius, not of `dx`.

The cell/Verlet unit suite uses exhaustive O(N²) reference graphs only for small cases, verifying force, density, conservative virial, viscous terms, capillary energy, and short trajectories. The native code also runs a one-time exhaustive audit for initial configurations smaller than 4,000 total particles. This audit is not used on large production runs.

## Validation scope

The standalone C++ regression suite passes, as does a C++ syntax-only check with local SPHinXsys/TBB interface stubs. These checks do **not** establish successful linking, correct VTP output, or performance on macOS. Only the on-device native CTest and A/B benchmark can confirm those claims. This change does **not** calibrate surface tension, contact angle, or capillary bridge force.
