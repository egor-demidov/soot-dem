# soot-dem + gas-free PF-SPH

The optional `pf_soot_dem` executable now uses the actual soot-dem `aggregate`
force functor, with libgran contact springs, bonded neck springs and Hamaker
forces. The old analytical coating-force functor is not added: fluid already
supplies the capillary load. The existing `pf_dem` executable remains available for its simpler model.
Fluid parameters are still provisional; calibration does not block these runs.

**Necks have no volume, mass, SPH boundary samples or wetting boundary.** They
only transmit DEM forces and moments. The fluid sees the moving primary-sphere
union. Optional neck VTP files contain diagnostic lines, not physical cylinders.

## Build against your soot-dem checkout

After merging this source archive into sphinxsys-capillary:

```bash
# Set this to the actual soot-dem repository location on your machine.
SOOT_REPO="$HOME/CLionProjects/soot-dem"

# Only needed if these dependency directories are missing:
git -C "$SOOT_REPO" submodule update --init deps/libgran deps/eigen

cmake -S . -B build-soot -DCMAKE_BUILD_TYPE=Release \
  -DPF_BUILD_SPHINXSYS_VALIDATION=OFF -DPF_WITH_SOOT_DEM=ON \
  -DSOOT_DEM_SOURCE_DIR="$SOOT_REPO"
cmake --build build-soot --target pf_soot_dem pf_soot_dem_test -j 6
ctest --test-dir build-soot -R '^pf_soot_dem_regression$' --output-on-failure
```

The build does not add soot-dem's entire CMake project, require its executables,
or build SPHinXsys. It compiles its force headers and Eigen into the coupled
executable. oneTBB is optional for the fluid pair loop; the mechanical loop is
serial. The optional historical SPHinXsys validation targets require an explicit
`PF_BUILD_SPHINXSYS_VALIDATION=ON` setting. They are off by default and do not
block configuration of this coupling. No bootstrap is called by CMake.

Tested source revisions: soot-dem `88ae15b513468aa352f0baf6784e9fb567a3bddf`,
libgran `ad9eb4d3214a2463aa555dae5293a041d96deb3b`, Eigen
`285da30ec309a0cacdb10d0fd15cc0de5dbc2f74`. A compatible local checkout can be
used; an unrecognized header API fails configuration rather than silently using
an invented substitute model.

## Run either coating initialization

```bash
python3 tools/run_pf_dem.py examples/soot_dem_droplet.xml \
  --executable "$PWD/build-soot/pf_soot_dem" --output results-soot-droplet

python3 tools/run_pf_dem.py examples/soot_dem_levelset.xml \
  --executable "$PWD/build-soot/pf_soot_dem" --output results-soot-levelset
```

These are development examples, not calibrated soot/water inputs. Both start
from the same aggregate and neck file. The levelset example has finer SPH spacing
so that its small separate menisci contain liquid samples. Coarser sampling can
leave this initialization empty; the runner reports that rather than creating
liquid elsewhere to compensate.

For only the short checks used during development:

```bash
python3 tools/check_pf_soot_dem.py --executable build-soot/pf_soot_dem \
  --output results-soot-smoke
```

The full example trajectories were not run locally. No equilibration, collision
benchmark or physical calibration is claimed by the short checks.

## Input

`<soot_dem>` replaces the generic `<contacts>` model. Do not enable both: that
would double-count dry contact forces. The XML runner disables generic contacts
automatically and rejects conflicting explicit input.

Its attributes use the existing soot-dem parameter names:

- Contact: `k_n`, `gamma_n`; `k_t`, `gamma_t`, `mu_t`, `phi_t`; corresponding
  `k_r/gamma_r/mu_r/phi_r` and `k_o/gamma_o/mu_o/phi_o`.
- Neck: `k_n_bond/gamma_n_bond`, `k_t_bond/gamma_t_bond`,
  `k_r_bond/gamma_r_bond`, `k_o_bond/gamma_o_bond`.
- `A`, `h0` for Hamaker force; `d_crit` for automatic INITIAL neck detection.
- `max_dt`: additional DEM timestep cap. The common adaptive step also respects
  stiffness/damping, fluid stability and rigid-motion limits.
- Optional `normal_yield`, `shear_yield`: uniform neck force thresholds using the
  existing soot-dem force-based breakage routine. Zero disables that criterion.
  Breakage is applied at the accepted endpoint; a trial/query cannot break a neck.

To set both thresholds from material strengths, give `neck_radius` in
meters and `normal_strength` and `shear_strength` in pascals. Supply all
three together as legacy XML `<let>` parameters or as modern `<soot_dem>`
attributes. The runner multiplies each strength by the neck area
`pi * neck_radius**2`. Use zero or omit the direct yield force inputs when
using strengths. The existing soot-dem routine breaks the neck when the
tensile elastic force or the sum of the sliding, rolling, and torsional
elastic force magnitudes reaches its threshold.

Initial topology is selected separately:

```xml
<necks file="soot_necks.csv"/>
<!-- or -->
<necks file="necks_0.vtk"/>
<!-- or upstream initial distance detection, including its cycle-prevention policy -->
<necks mode="auto"/>
<!-- or no initial bonds -->
<necks mode="none"/>
```

CSV uses zero-based `i,j` columns. Legacy soot-dem ASCII neck VTK files are read
from their encoded `connections` field. Explicit graphs may contain cycles;
automatic graph creation retains upstream's forest policy. Broken necks do not
regrow during dynamics. Spatial contact can still form after breakage, through
the nonbonded soot-dem contact law.

The current soot-dem aggregate functor assumes equal primary radii and masses.
The adapter enforces that restriction; it does not silently apply monodisperse
coefficients to a polydisperse aggregate. Particle fixity is supported.

Particle inputs remain CSV or ASCII VTK/VTP. For a legacy soot-dem particle VTK
whose coordinates are normalized by primary radius, an SI example is:

```xml
<particles file="aggregate.vtk" radius="1" density="1800"
           length_scale="14e-9" ignore_velocities="true"/>
```

Both the fallback radius 1 and positions are multiplied by length_scale. Supply
fluid, liquid geometry and DEM coefficients in the same consistent unit system;
the example model-unit coefficients are not SI material values.

Coating geometry remains:

```xml
<liquid mode="droplets">
  <droplet x="0" y="0" z="0" radius="0.6"/>
</liquid>
<!-- alternatively -->
<liquid mode="levelset" alpha="0.25" threshold="0.6" rm="0.24"/>
```

Droplets fill their union minus primary spheres. The levelset fills free-space
samples where sum_i 0.5 erfc((distance_i/rm - 1)/alpha) reaches the threshold.
rm=0 or omitted uses each sphere radius. Necks do not participate in either
calculation. After initialization the liquid evolves freely; the levelset is
not reapplied as particles move and does not add or remove liquid.

## Coupling and history integration

Native soot-dem pair responses are linear/angular ACCELERATIONS. The adapter
converts them to forces and moments using the dry particle mass and inertia.
It evaluates one canonical pair and assembles its reciprocal force/moment at the
other endpoint. This agrees with the native pair response for equal spheres and
preserves total pair angular momentum without a numerical projection.

Each accepted step uses:

1. A nonmutating initial force query from committed DEM history.
2. A temporary half-step spring-history predictor from the old state.
3. A midpoint fluid/DEM force evaluation from that trial history.
4. A full history update from the OLD committed state using midpoint kinematics,
   followed by endpoint neck breakage and commit.

The full step does not add a second half-step onto the trial predictor. Queries,
CSV output, rejected steps and timestep estimation cannot age the history.
Particle positions, quaternions and liquid states use the same synchronized
integrator; the old soot-dem integrator is not also advanced independently.

## Small upstream overlay

`tools/prepare_soot_backend.py` copies only the required headers into
`build-soot/soot-backend`. It never edits the source checkout. The overlay adds
spring-state access and variable spring-update dt for trial/commit and restart.
It also removes an unused granular-integrator include from aggregate.h.

One constitutive implementation correction is included and tested: torsional
spring history must be projected onto the pair axis, while sliding/rolling
histories remain tangent to it. The existing common projection erased accumulated
axial torsion on the next force evaluation. The overlay fixes that projection in
both contact and neck functors. Other upstream constitutive behavior is retained.
`backend.json` records source revision, content fingerprint and modifications.

## Restart and output

```bash
python3 tools/run_pf_dem.py examples/soot_dem_droplet.xml \
  --executable "$PWD/build-soot/pf_soot_dem" \
  --resume "$PWD/results-soot-droplet/dem_checkpoint.dat" \
  --end-time 1 --output results-soot-droplet-extended
```

PFDEM version 2 checkpoints include both fluid/body state and the mechanical
backend: current graph including broken necks, bonded and nonbonded spring
histories, parameters and source fingerprint. Earlier PFDEM version 1 checkpoints
without a mechanical backend remain readable. Restart rejects a changed physical
XML configuration or changed soot-dem backend code. It does not regenerate the
neck graph or heal broken necks. If backend code changes, retain the old build for
exact continuation or start a new run explicitly.

Existing fluid/sphere outputs and surface exporter remain usable. New output:

- `dem_particles.csv` has `dem_fx/fy/fz/tx/ty/tz` separately from fluid loads.
- `dem_history.csv` includes `neck_count`.
- `dem_neck_forces.csv` records one row per active neck at each DEM CSV
  recording step, including the initial and final records. It contains
  `step,time,particle_i,particle_j,normal_force_N,sliding_force_N,rolling_force_N,torsional_force_N`.
  Time is in seconds and all four spring forces are in newtons. Normal
  force is signed, positive in tension. Sliding, rolling, and torsional
  forces are the magnitudes used in the breakage criterion, each equal
  to its bond stiffness times its stored spring magnitude. During a
  fluid hold, they report the retained elastic state. Broken necks
  disappear from subsequent records. The pair of particle IDs identifies
  a neck across frames and checkpoint restarts.
- `necks_*.vtp` displays current bonds as lines when dumps are enabled.

`builtin_contact_energy` reports only the generic model's Hertz energy and is
zero for this backend. Total soot-dem elastic/Hamaker energy is NOT yet exported;
momentum conservation does not substitute for an energy-convergence test.

## Verification

Short C++ checks compare adapter force/moment conversion against the actual
soot-dem functor, test invisible neck geometry, repeated queries and rollback,
exact one-dt sliding-history accumulation, persistent axial torsion, fixed-support
momentum, tensile breakage and non-healing restart. Coupled trajectories compare
exactly across a split/restart, including spring histories and quaternions.
Short executable checks cover both coating modes and their restarts. Python
checks cover topology imports, contact double-count prevention and uniform-radius
validation. No long physical run was performed.

C++20 executables were compiled directly with g++ against the real downloaded
headers. CMake was unavailable locally, so configuration should be checked on
your machine using the commands above. Dense upstream N-by-N histories and serial
DEM pair evaluation currently limit scalability. Sub-SPH gap lubrication and
quantitative wetting/force validation remain unresolved as before.
