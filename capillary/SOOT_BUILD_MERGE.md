# Add the coupling to the original soot-dem build

The standalone `CMakeLists.txt` in earlier archives was not a replacement for
soot-dem's root build. This installer preserves the original soot-dem target list
and adds the fluid coupling as a `capillary/` component.

## Install / repair

Extract the archive into a temporary directory, then point the installer at the
ACTUAL soot-dem source repository (the directory containing `src/aggregate.h`
and `src/restructuring.cpp`):

```bash
unzip -o ~/Downloads/sphinxsys_capillary_neighbor_acceleration_full_source.zip \
  -d /tmp/sphinxsys-update

python3 /tmp/sphinxsys-update/sphinxsys-capillary/tools/install_into_soot_dem.py \
  --soot-repo "$HOME/CLionProjects/soot-dem"
```

Change the repository path if your checkout is elsewhere. Do not copy this
archive's standalone root CMakeLists over soot-dem's root again.

The installer:

- Preserves an existing original root build, including custom targets/options.
- If that file was overwritten by the standalone capillary build, recovers the
  most recent original soot-dem root build from local Git HEAD/history.
- Backs up the current root file as `CMakeLists.txt.before-capillary-merge`
  before changing it, retaining numbered backups when needed.
- Places only the declared component files in `soot-dem/capillary/`.
- Adds `cmake/SootCapillaryCoupling.cmake` and one include at the end of the root
  build. It does not remove original sources, dependencies, builds or runs.

If no original build can be found in local Git history, installation stops before
replacing the root file. `cmake/soot-dem.CMakeLists.merged` is a complete merged
REFERENCE based on upstream soot-dem commit 88ae15b; it retains all 24 original
executables. The installer preferentially uses your own original file so custom
changes are preserved rather than replaced by that reference.

## Configure and build the ORIGINAL repository

```bash
cd "$HOME/CLionProjects/soot-dem"
cmake -S . -B build -DPF_BUILD_SPHINXSYS_VALIDATION=OFF
cmake --build build --target pf_soot_dem pf_soot_dem_test -j 6
ctest --test-dir build -R '^pf_soot_dem_regression$' --output-on-failure
```

Use your usual CLion/CMake build directory and compiler settings if different.
If an old build directory belongs to a different SOURCE directory, select a new
build directory rather than reusing that incompatible CMake cache.

The original restructuring, aggregation, AFM, mechanics, validation and neck
parametrization targets remain available. The new targets are added alongside
those. `SOOT_BUILD_CAPILLARY_COUPLING=OFF` disables just the added component.
SPHinXsys validation remains an optional explicit opt-in, not a bootstrap
requirement. The containing soot-dem repository supplies libgran and Eigen;
there is no need to set a separate SOOT_DEM_SOURCE_DIR for the merged build.

## Run from the original repository root

```bash
python3 capillary/tools/run_pf_dem.py capillary/examples/soot_dem_droplet.xml \
  --executable "$PWD/build/pf_soot_dem" --output results-soot-droplet

python3 capillary/tools/run_pf_dem.py capillary/examples/soot_dem_levelset.xml \
  --executable "$PWD/build/pf_soot_dem" --output results-soot-levelset
```

The model/input details remain in `capillary/SOOT_DEM_COUPLING.md`. Input-relative
aggregate and neck paths continue to work. This build integration does not alter
the force laws, spring histories, checkpoint format or fluid initialization.

## Checked in this update

Configured the actual upstream soot-dem project with the component installed.
All 24 original executable targets and the added coupling targets were available.
Built the coupled backend regression through CMake and ran it through CTest: pass.
The original targets were preserved, not rebuilt or benchmarked. Installer checks
cover custom target preservation, recovery from Git after an accidental overwrite,
and repeated installation. No simulation trajectory was run for this build fix.
