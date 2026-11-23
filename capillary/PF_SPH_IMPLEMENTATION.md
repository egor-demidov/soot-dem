# Gas-free PF-SPH: independent pair module

Status (2026-11-03): the force module and liquid-only adapter are implemented.
See **PF_LIQUID_RUNS.md** for the new `pf_capillary_validation` executable,
operator details, validation limitations and the two-run pilot script.
The existing `capillary_validation` executable remains two-phase.

## Published force law selected for this implementation

Whebell, Moroney, Turner, Pethiyagoda & McCue, JFM 1019, A37 (2025),
https://doi.org/10.1017/jfm.2025.10601, equations (2.19)-(2.21).
The paper supplies a liquid-only formulation with separate fluid-fluid and
fluid-solid polynomial profiles. This is not the Gaussian law of Tartakovsky &
Panchenko (2016), nor the cubic-spline law used by Shigorina et al. (2017).
We do not combine their coefficients or import their surface-tension mappings.

For r_ij=x_i-x_j, r=|r_ij|, q=r/H and common support H:

    F_ij = H s_ij f_ij(q) r_ij/r,
    f(q) = (1-q)^2 [1+2q+(3+alpha)q^2]       for 0 <= q < 1,
           0                               for q >= 1.

The published alpha values are -23.5 for fluid-fluid and -11 for fluid-solid.
With the paper's H=4*dx, the repulsion/attraction zero crossings are approximately
1.1*dx and exactly 2*dx. Positive f is repulsive. The implementation intentionally
has no contact-angle input: s_ff and s_fs are independent nonnegative strengths.
Their units are N/m in SI, but s_ff is NOT automatically the measured surface tension.
Uniform sampling and fixed support are assumed; variable-resolution quadrature is
not included or validated.

## Shared potential and units

The pair energy is derived by integrating the published force. Setting U(H)=0,
with z=1-r/H gives

    U(r) = s H^2 z^3 [(alpha+6)/3 - (alpha+4)z/2 + (alpha+3)z^2/5], r < H,
           0, r >= H.

Then -dU/dr = H*s*f(r/H). The code evaluates this factored expression rather than
subtracting two nearly equal antiderivatives near the cutoff. Energy has units J;
force has units N. Division by a particle mass occurs only when the caller forms
acceleration, not when assigning action and reaction.

Self pairs are excluded by neighbor traversal. Coincident distinct active samples
raise an explicit error because the published radial law has no direction there.
The module does not invent a direction, cap a force, or hide the overlap. This does
not forbid DEM particle overlap; samples inside the overlapping solid union require
an exposed-surface strategy in the future adapter.

## Local reciprocal transfer and rigid loads

`src/pf_capillary.h` is independent of SPHinXsys and all third-party headers.
`evaluate` returns the force on i, pair energy and scalar pair virial. A unique
undirected pair is evaluated once; `accumulate` adds the same force to i and its
negative to j. Pair type, support and strength must be symmetric under i/j exchange.
The caller must not include a pair twice or apply this capillary model to solid-solid
DEM contacts. Those contacts remain the DEM model's responsibility.

For a boundary sample j attached to body k, store its force before reducing:

    F_k = sum_j F_j,
    T_k = sum_j (x_j-X_k) cross F_j.

Central reciprocal pair forces give both zero total force and zero total torque
about a common origin, up to floating-point accumulation. This includes the fluid
moment and the body's orbital and spin moments. It does not mean a body's individual
torque is zero, nor does it prove every other SPH operator has those identities.
A homogeneous sphere's continuum conservative capillary torque has additional
symmetry constraints; discrete sampling must not be mistaken for physical roughness.

Finite timestep energy conservation is not claimed. Even exact potential gradients
have time integration error. Viscosity dissipates energy and needs a stable update;
retain the acoustic-viscosity correction established in the current adapter.

## What passed locally

The standalone C++ test remains active in optimized Release/NDEBUG builds and checks:

- Published support, signs, zero crossings and smooth cutoff.
- Analytic potential derivative against the actual force.
- Full multi-sample energy gradient by coordinate finite differences.
- Global force and torque reciprocity without a residual correction.
- Energy invariance and force covariance under rigid frame transformations.
- Rigid translation and rotation virtual work with a nonzero body torque.
- Unequal-mass pair impulses, and dimensional scaling of force and energy.
- Invalid support and coincident active sample rejection.

These are operator checks, NOT surface-tension, contact-angle, bridge-force or
fluid-solver validation. Test strengths are arbitrary and are not production inputs.

Quick local build, requiring only a C++17 compiler:

```bash
c++ -std=c++17 -O2 -DNDEBUG tests/test_pf_capillary.cpp -o /tmp/pf_capillary_test
/tmp/pf_capillary_test
```

Or configure only the independent module, without SPHinXsys dependencies:

```bash
cmake -S . -B build-pf-pairs -DCAPILLARY_PAIR_TESTS_ONLY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-pf-pairs --target pf_capillary_test
ctest --test-dir build-pf-pairs --output-on-failure
```

## Integration safeguards (original pair-module plan)

Items 1–4 now have an initial implementation in `pf_validation.cpp`; see
PF_LIQUID_RUNS.md for what was actually tested. Physical calibration and moving
DEM validation remain outstanding.


1. Construct a genuinely liquid-only adapter using SPHinXsys storage/neighbors/I/O;
   pressure, continuity, viscosity and solid boundary treatment must form a compatible
   scheme. Do not just remove air from the old multiphase operator sequence.
2. Use support H=4*dx consistently for this selected profile. The current library's
   default neighbor radius is smaller; using it would silently truncate interactions.
   Rebuild/validate neighbor lists often enough that no entering pair is omitted.
   Force values must use current positions, not stale cached distances.
3. Keep all reciprocal boundary contributions. Pressure and viscous forces need the
   same audit as capillarity. Do not infer solid loads from unrelated SPH acceleration
   reductions or bolt on a guessed contact-line load.
4. Disable old CSF and prescribed-normal wetting when enabling this model. Maintain
   separate archived baselines rather than adding the models together.
5. Calibrate actual mechanical surface tension. Pair forces contribute to mechanical
   stress; an EOS-pressure-only Laplace comparison can be wrong. The scalar virial
   returned by the module is an ingredient, not a ready pressure estimator: bulk
   averaging, domain volume, EOS/kinetic terms and boundary selection still matter.
6. Check planar/spherical interfaces, a flat-wall droplet and then a two-sphere bridge.
   The literature's strength-to-angle/surface-tension relationships cannot be assumed
   exact for different kernels, EOS, resolution or solid boundary sampling.
7. Only then enable moving DEM, check impulse/torque/energy accounting, and test
   bridge joining/breakup. Gas pressure/cushioning and trapped-bubble compression
   remain omitted physical effects.

The pair module remains portable; the separate liquid-only adapter is an
unvalidated solver baseline and does not claim moving-DEM validation for the
2025 model. Existing triangle virtual-work tests remain independent reference
checks; no runtime interface mesh has been introduced.
