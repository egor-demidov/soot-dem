# General moving-particle DEM / PF-SPH coupling

The actual soot-dem backend is now available separately: [SOOT_DEM_COUPLING.md](SOOT_DEM_COUPLING.md).

This stage enables development with provisional material parameters. It does not
wait for the surface-tension/contact-angle calibration. Calibration reports retain
their actual failed status; none of the inputs is relabeled as validated water.

## What is implemented

`pf_dem` runs an arbitrary assembly of spherical DEM particles coupled to the
same gas-free PF-WCSPH equations used by `pf_capillary_validation`. It reuses the
existing capillary pair law, pressure, continuity, viscosity, deterministic pair
ordering and cached Verlet search. It does not add a second fluid model.

The general executable is independent of SPHinXsys storage/I/O, so it can also
be built without that dependency. Existing SPHinXsys validation cases are retained
unchanged. If oneTBB is discoverable, `pf_dem` uses the existing deterministic
parallel force evaluator; otherwise it explicitly reports a serial fallback.

- Every primary particle has its own dry mass, center, velocity, laboratory angular
  momentum, and scalar-first body-to-world quaternion.
- Boundary samples have immutable local coordinates. At every force stage they
  follow the body's translation and rotation, with velocity V + omega cross r.
  Their quadrature masses are NOT added to the particle's physical mass/inertia.
- Each fluid pair is evaluated once. The solid endpoint's exact opposite force
  is reduced to its owner; torque is the cross product of that sample's lever arm
  with its force. There is no separate contact-line estimate or global projection.
- All bodies and liquid advance on a synchronized explicit midpoint step. Final
  positions use trapezoidal velocity. The contact, acceleration, sound-speed,
  viscous and rotation bounds constrain the common timestep.
- Contacts allow overlap: Hertz elasticity, optional normal dashpot, and optional
  regularized kinetic Coulomb friction. Tangential reactions act at a common point
  inside the overlap, so their moments preserve pair angular momentum.
- Boundary samples buried inside another particle are excluded from fluid pairs.
  IDs and local coordinates remain unchanged, so there is no cell remapping,
  fluid deletion, donor search, or refill operation when solids move.
- Fixed particles remain fixed; their impulses/moments go into a support ledger.
- The `PFDEM 2` checkpoint (with legacy PFDEM 1 reading) stores liquid states, body states, sample ownership and
  local coordinates, clock, contact parameters and momentum ledgers. Resume
  reconstructs the neighbor cache and checks the prepared input signature.

The built-in contact model is a working DEM model, but is NOT the soot-dem bonded
neck model. An interface for supplying that model is provided below. No libgran
or libtimestep linkage is claimed by this package.

## Build and short checks

An existing configured native build can build the added targets:

```bash
cmake --build build --target pf_dem pf_dem_test -j 6
ctest --test-dir build -R '^pf_dem_regression$' --output-on-failure
python3 -m unittest discover -s tests -p test_pf_dem_runner.py -v
python3 tools/check_pf_dem.py --executable build/pf_dem --output results-dem-smoke
```

Alternatively, build this driver without SPHinXsys:

```bash
cmake -S . -B build-dem -DCMAKE_BUILD_TYPE=Release -DCAPILLARY_PAIR_TESTS_ONLY=ON
cmake --build build-dem --target pf_dem pf_dem_test -j 6
ctest --test-dir build-dem -R '^pf_dem_regression$' --output-on-failure
python3 tools/check_pf_dem.py --executable build-dem/pf_dem --output results-dem-smoke
```

The check script runs only tiny startup/restart trajectories. It does not equilibrate
an interface, complete a collision, or validate a capillary bridge force.

## Run the two-particle example

```bash
python3 tools/run_pf_dem.py examples/dem_two_particles.xml \
  --executable build/pf_dem --output results-dem-two
```

One sphere initially intersects the droplet. The second starts outside both the
liquid and the finite force-support range, moving along +x toward it. Its initial
fluid load is exactly zero. The full example was NOT run locally in this stage;
it remains a development run. Use fresh output directories for each run.

To inspect only initialization, without starting dynamics:

```bash
python3 tools/run_pf_dem.py examples/dem_two_particles.xml \
  --output prepared-dem-two --prepare-only
```

## Aggregate / levelset example

```bash
python3 tools/run_pf_dem.py examples/dem_aggregate_levelset.xml \
  --executable build/pf_dem --output results-dem-levelset
```

Edit `<particles file="...">` to load a CSV (`x,y,z` plus `radius`, or a default
radius attribute), legacy ASCII VTK POINTS, or single-piece ASCII VTP Points.
Optional point radius arrays are named `radius` or `radii`. These are particle
centers, not a triangulated surface mesh. Binary VTK/VTP is rejected with an
explanation. `length_scale` scales imported centers/radii. Density uses the fluid
model's consistent units. File velocities are ignored by default, as requested;
CSV velocities may be enabled explicitly with `ignore_velocities="false"`.
Inline `<particle>` elements may instead prescribe vx/vy/vz and wx/wy/wz.

Liquid initialization supports:

- `<liquid mode="droplets">` containing one or more spherical `<droplet>` elements:
  fill their union after subtracting the union of all solid particles.
- `<liquid mode="levelset" alpha="0.3" threshold="0.5" rm="0.25"/>`:
  fill free-space samples for which sum_i 0.5 erfc((distance_i/rm - 1)/alpha)
  is at least the threshold. Omit rm, or set rm=0, to use each particle's radius.

Bounds enclosing particles and liquid plus `domain margin` are computed
automatically. For levelsets, a conservative erfc-tail bound is used. The bounds
only set the INITIAL sampling region; they are not periodic walls and do not clip
subsequent motion. `max_samples` also guards excessive initialization work.
There is no explicit gas phase and no evaporation or condensation source.

The example aggregate is three overlapping, UNBONDED spheres. It will respond to
Hertz contact and fluid forces; it is not a permanently connected soot backbone.

## Restart

Keep the same XML physical settings and aggregate file; runtime/output frequency
changes are allowed. The target is an absolute simulation time. `--steps` limits
ADDITIONAL steps of this invocation. A fresh output directory is required:

```bash
python3 tools/run_pf_dem.py examples/dem_two_particles.xml \
  --executable build/pf_dem --resume results-dem-two/dem_checkpoint.dat \
  --end-time 4 --output results-dem-two-resumed
```

Restart rejects changed fluid parameters, geometry, velocities, body material,
fixity or contact settings. It never silently resets fluid density or regenerates
liquid particles. Histories in a resumed folder start at the checkpoint time.

## Output

Frequencies are XML run attributes: `record_every_steps`, `dump_every_steps`,
`checkpoint_every_steps`. Setting dump/checkpoint frequency to zero disables
periodic writes; a final checkpoint is still written. `threads` is honored when
TBB is available. No OpenMP implementation is added to this driver.

- `dem_particles.csv`: body motion, quaternion, and separate fluid/contact loads.
- `dem_history.csv`: constant particle mass, compressible liquid volume sum(m/rho),
  momentum errors including supports, kinetic/contact energy and neighbor counts.
- `liquid_*.vtp`: liquid point clouds with density, pressure and velocity. These
  are NOT exported isosurfaces; use the existing point-cloud reconstruction.
- `particles_*.vtp`: triangulated body surfaces. Color by `ParticleId` to identify
  particles, or `MaterialZ` to see rotation of the material surface.
- `dem_checkpoint.dat`, `resolved_input.json`, `prepared.dat`, `run.log`.

All examples use one self-consistent set of MODEL units. The runner does not
convert a calibrated macroscopic gamma/contact angle into sff/sfs. For SI inputs,
all lengths, densities, speeds, viscosity and contact coefficients must use SI
consistently; provisional strengths still require calibration at the chosen dx.

## Integrate an external DEM law

`pf_dem::System` exposes `bodies`, `sites`, and `Stage::hydrodynamic`. To use the
soot-dem contact/neck law, supply a pure force evaluation callback:

```cpp
system.contact_model.stiffness = 0.; // disable built-in contacts when replacing
system.contact_model.damping = 0.;
system.contact_model.friction = 0.;
system.external_max_dt = dem_stable_step;
system.extra_loads = [&](const std::vector<pf_rigid::Sphere>& stage_bodies,
                         double stage_time) {
    // Convert synchronized centers, velocities, omega and quaternion to DEM.
    // Return ONE force/torque per body, from contacts/bonds/external fields only.
    // Never return fluid loads here: those are already added by the coupler.
    return evaluate_dem_laws(stage_bodies, stage_time);
};
```

The coupler performs the integration; do not advance a second DEM integrator in
the callback. Stateful tangential springs/breakage need trial/commit handling and
their own checkpoint serialization. Saving with a callback currently throws
rather than silently discarding that external law's state. The callback's net
applied impulse is included in conservation diagnostics. Built-in Hertz/friction
restart is fully implemented; use the stateful backend in SOOT_DEM_COUPLING.md for actual soot-dem coupling.

## Verification and practical limits

Local checks: g++ C++17 compilation; N-body regression with reciprocal forces and
friction torques; fixed supports; exact split trajectory restart; callback stages;
buried sample masking; solid-union liquid subtraction; four XML/import checks;
and tiny executable startup/restart checks. See `verification/dem_coupling/`.
No long collision or equilibration run, TBB build, or physical calibration was
performed locally. CMake was unavailable locally; the driver and core tests were
compiled directly with g++.

Dry contacts and body-overlap screening currently use O(N^2) body-pair loops;
fluid interactions use the cached spatial search. Occlusion tests visit only
neighbors whose solid spheres overlap. Very thin films/gaps below the SPH spacing
and large overlaps are not physically resolved. The binary masking of buried
quadrature may change forces abruptly near contact; mass and momentum checks do
not establish energy convergence through those topology changes. Compressibility
also means constant liquid particle mass does not imply exactly constant volume.
Bonded necks, calibrated dynamic wetting, lubrication below dx, and quantitative
force/energy convergence remain separate work.

## Export a smooth liquid surface from an existing run

The offline exporter uses the same mass-weighted Wendland reconstruction as
`pf_interface.py`, contouring reconstructed density at rho0/2. It does not contour
the nearly uniform density values attached directly to the particle points.
It writes triangulated VTP files with smooth normals and a ParaView PVD collection.

```bash
python3 -m pip install -r requirements.txt
# Preview just the final saved frame first:
python3 tools/export_dem_surface.py results-dem-two-2 \
  --last-only --output results-dem-two-2/surface-preview
# Then export the saved trajectory, or add --every 5 to reduce export work:
python3 tools/export_dem_surface.py results-dem-two-2 \
  --grid-dx 0.5 --output results-dem-two-2/liquid-surfaces
```

Open `liquid-surfaces/liquid_surface.pvd` in ParaView and display as Surface.
Hide the liquid point cloud and keep the particle surfaces opaque. Normals are
included for smooth shading. A light liquid color and moderate specular lighting
can improve contrast. For closer views, `--grid-dx 0.25` gives a finer VISUAL grid;
it does not increase the simulation resolution. Keep a fixed grid setting when
comparing frames. Existing export directories are preserved; choose fresh outputs.

The reconstruction can extend inside a solid near contact. Solid geometry is not
clipped into this exporter: opaque particle surfaces hide that portion. Interfaces
and thin menisci near solids remain sensitive to SPH spacing and support. Surface
mesh enclosed volume is not a conservation diagnostic; inspect `liquid_mass` and
`liquid_volume` in the simulation history separately.

Most default dump steps have exact times in `dem_history.csv`. If selected frames
lack times there or in VTP FieldData, the PVD explicitly uses step numbers for the
whole sequence rather than inventing simulation times. A resumed run is exported
from its own folder with its retained absolute times.

The exporter needs NumPy, SciPy and scikit-image (already in requirements.txt).
Local checks covered compact-kernel mass integration, VTP I/O, normals, frame
selection and series metadata. scikit-image was unavailable locally, so actual
marching-cubes extraction was not executed in this update. No simulation was run.
The relative-executable-path bug in the XML runner is also fixed and regression-tested.
