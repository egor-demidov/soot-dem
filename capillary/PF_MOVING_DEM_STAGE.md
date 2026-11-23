# Stage: moving rigid sphere with reciprocal liquid-only PF-SPH loads

**Status:** Core 6-DOF midpoint mechanics and unit checks are implemented.
The `pf_capillary_validation --case immersed` SPHinXsys adapter is **experimental**:
native SPHinXsys/TBB linking and short trajectories must be run against your pinned
checkout on your machine. The sandbox lacks the native checkout and dependencies.
Do **not** interpret these tests as validation of a wetting contact angle or capillary bridge.

## Implemented

- `src/pf_rigid.h`: freely translating/rotating rigid spheres with physical dry
  mass `rho_s*(4*pi/3)*R^3`, isotropic spherical inertia `2mR^2/5`, body-to-world
  unit quaternion, *world-frame angular momentum* and per-body force/torque.
- Every PF fluid–solid pair is evaluated once in `pf_liquid::accumulate_pair`.
  Fluid force is `F`, solid sample force is `-F`. The resulting rigid load
  is `sum_a f_a`, torque `sum_a [(y_a-X) cross f_a]` for samples belonging to
  a body. Solid quadrature-point masses are **not** added to the dry DEM mass.
- At midpoint and end of the step, all solid boundary samples have their
  positions and velocities regenerated as `X + R(q)r_a` and
  `V + omega cross R(q)r_a`, with `omega=L/I`. The native fluid and moving
  solid neighbor cell-linked lists are refreshed at evaluation positions.
- The experimental native `immersed` case places a moving sphere inside a
  finite spherical volume of liquid without gas and with no gravity. The
  initial fluid lattice excludes the sphere core. The sphere begins with
  optionally nonzero x-velocity and z-spin to exercise real motion.
- Standalone sphere–sphere contact: non-dissipative normal force
  `k * overlap^(3/2)` and potential `(2/5) k * overlap^(5/2)`.
  **This is NOT yet enabled in the single-sphere native scenario.**
- The full-state checkpoint is upgraded to `PFSPH_CHECKPOINT 2`, storing sphere
  positions, velocities, angular momenta, orientations, material parameters,
  initial conserved total linear/angular momentum, wall impulse ledger and
  all fluid/solid samples. Old v1 files remain readable. When restarting
  a moving sphere, its radius, mass, sample order and transformed boundary
  geometry must match the initialized geometry.
- `pf_rigid.csv` now reports rigid kinematics, SPH force and SPH torque;
  `pf_history.csv` linear/angular errors for `immersed` compare TOTAL
  liquid+DEM momentum against conserved initial values.
- `tools/run_pf_moving_smoke.py` runs a very short native adapter check,
  checks finite diagnostics and sphere motion, and packages the log, history,
  checkpoint and performance measurements in a ZIP.

## Compile and run the *short* native check

Update the files over the existing project source. Keep your existing
`deps/SPHinXsys` checkout (do not replace it). From the project root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target pf_capillary_validation pf_rigid_test pf_restart_test -j 6
ctest --test-dir build -R 'pf_(capillary_pair|liquid|restart|snapshot|rigid)_regression' --output-on-failure
ctest --test-dir build -R '^pf_immersed_adapter_smoke$' --output-on-failure
```

Then collect a clean, isolated 0.002-unit run:

```bash
python3 tools/run_pf_moving_smoke.py \
  --executable "$PWD/build/pf_capillary_validation" \
  --output results-pf-moving-smoke
```

Only send `results-pf-moving-smoke.zip` if the CTest smoke succeeded;
otherwise send the complete native compiler/runtime error first.

Changing `--radius` changes the outer liquid radius; `--sphere-radius`
changes the solid sphere, with `outer_radius > sphere_radius + 2dx`.
`--sphere-density` controls the *physical* DEM mass and inertia (not the
volumetric quadrature weights of PF boundary samples).

To restart a moving-sphere run from its **v2 checkpoint** (always in a new
work directory):

```bash
../build/pf_capillary_validation --case immersed \
  --dx 0.12 --radius 0.8 --sphere-radius 0.24 --sphere-density 2 \
  --sff 0.001 --sfs 0.002 --viscosity 0.05 --sound-speed 10 \
  --end-time 0.004 --dump-every 0 --record-every 1 \
  --checkpoint-interval 0.002 \
  --resume ../results-pf-moving-smoke/pf_checkpoint.dat
```

The initial `--sphere-vx`/`--sphere-omega-z` are ignored when loading a
full checkpoint, which replaces sphere positions, momenta and orientation.

## What is *not* implemented or validated

1. No two-sphere *native* bridge, no free-interface convergence/meniscus
   measurement, and no physical surface tension/Young contact-angle calibration.
   `sff`, `sfs` are raw PF strengths and do not equal gamma.
2. No calibrated solid impermeability/contact boundary treatment. At long
   times, particles may penetrate or the PF model may be unstable; test
   geometric overlap and density, especially with moving solid boundaries.
3. No history-dependent DEM tangential friction, contact damping, bonded
   bridge transfer, or rupture criterion. Native `immersed` has one sphere;
   only a standalone Hertz-type two-sphere normal-contact law is tested.
4. No interface to your existing soot DEM C++ library yet. The new core is
   deliberately free of SPHinXsys, making it suitable for a future adapter.
5. Native execution and scaling, capillary bridge force/torque vs Young–Laplace
   and refinement tests remain **unverified**. The independent source tests
   check discrete reciprocity, virtual-work capillary force/torque derivatives,
   finite-state updates and exact split/restart in a fixed-order graph.

## Results checked here

- Standalone: 5 C++ CTests passed; 22 Python unit tests passed.
- Native adapter: **syntax-only** compiled against minimal temporary
  SPHinXsys/TBB interface stubs, NOT linked against the real upstream code.
- A true SPHinXsys build and runtime smoke must be completed on the user's
  host before attempting any long bridge or wetting simulation.
