# Capillary force and torque: consistent closure and implementation gates

Status, 2026-11-03: derivation and independent geometric reference implemented;
production SPH capillarity and wall loads are unchanged. The successful acoustic
viscosity update remains enabled. No new fluid benchmark is requested at this stage.

## 1. What the measured defect establishes

The t=2 acoustic-viscosity audit reports a supported momentum defect of 0.051140,
and a maximum residual of 0.00006783 when the applied continuum surface-force (CSF)
impulse is accounted for as an external load. This identifies the dominant term in
the budget, but does not determine a unique physical force distribution on solids.
It does not show that every CSF contribution belongs to a nearby particle.
Pressure wall extrapolation already depends on prior acceleration, including CSF;
adding a separately guessed capillary force can double-count parts of the response.

The pinned upstream wetting example rotates liquid normals and applies a curvature
force; it does not differentiate a common fluid/solid interfacial energy. Its liquid
and gas pressure operators also differ, and solid pressure loads are sampled at a
different time level from fluid updates. A capillary-only improvement will not make
all of those other operators conservative automatically.

## 2. Continuum energy and Young's condition

For an isothermal, homogeneous, rigid solid with constant total exposed surface
area, the position-dependent interfacial energy, up to a constant, is

    E_cap = gamma_LG A_LG + (gamma_SL - gamma_SG) A_SL
          = gamma [A_LG - cos(theta_Y) A_SL].

The contact angle theta_Y is measured through the liquid. The second equality is
Young's relation, gamma_SL-gamma_SG = -gamma cos(theta_Y). This representation
requires neither prescribed advancing/receding angles nor a hysteresis model.
Surface chemistry, hysteresis and dissipation would require additional terms or
constitutive laws; they must not be inferred from a discretization-dependent angle.

Consider a flat wall, and let e_r point outward from the wetted patch and e_z into
the fluid. The outward co-normal of the liquid-gas surface at its circular edge is
m = cos(theta) e_r - sin(theta) e_z. The boundary force associated with liquid-gas
area variation is -gamma*m. Wetting-area variation adds
-(gamma_SL-gamma_SG)*e_r. Their combined tangential coefficient is

    -[gamma cos(theta) + gamma_SL-gamma_SG].

It vanishes at theta=theta_Y. The normal component on the solid is upward,
gamma sin(theta) per unit contact-line length. These signs follow from virtual
work; they are not selected by fitting a measured global momentum residual.

For a spherical cap of curvature radius R and contact radius a=R sin(theta),

    F_line,z = 2 pi a gamma sin(theta),
    F_pressure,z = -(2 gamma/R) pi a^2,
    F_solid,z = 0.

Thus a fixed weightless equilibrium droplet can have a substantial downward pressure
load and equal upward contact-line load while its total load on the plate is zero.
This explains the scale and sign of the observed missing budget term without
establishing that the existing discrete CSF is exactly the continuum curvature load.
A bridge, curved wall or accelerating droplet generally does NOT have zero net load.

## 3. Discrete energy, reciprocity and rigid bodies

Let x_i denote fluid positions and let solid samples be y_ka=X_k+R_k r_ka.
A capillary discretization must define one differentiable energy

    E_h(x,y) = gamma A_LG,h(x,y) - gamma cos(theta_Y) A_SL,h(x,y).

Forces are defined by differentiating this same expression:

    f_i = -dE_h/dx_i,       g_ka = -dE_h/dy_ka,
    F_k = sum_a g_ka,
    T_k = sum_a (y_ka-X_k) cross g_ka.

Equivalently, rigid translations and infinitesimal rotations give

    delta E_h = -sum_i f_i dot delta x_i
                -sum_k [F_k dot delta X_k + T_k dot delta phi_k].

If the energy is invariant under simultaneous translation of every sample,

    sum_i f_i + sum_k F_k = 0.

If it is also invariant under simultaneous rotation,

    sum_i x_i cross f_i + sum_k [X_k cross F_k + T_k] = 0.

These are algebraic force identities before time integration, not corrections.
An antisymmetric pair force alone guarantees only linear momentum. Angular
momentum additionally requires central pair forces, or a complete rotationally
invariant multi-body energy derivative. A general symmetric surface-stress tensor
in the continuum does not by itself prove its SPH discretization has this property.

The same force and torque must enter the rigid-body integrator and the balance
ledger at the same force evaluation. Comparing moments about different origins or
sampling lever arms at different times can create an apparent defect. No scalar
choice of a DEM integrator repairs a non-reciprocal force operator.

The identities do not imply exact finite-step energy conservation. For example,
symplectic Euler has integration error even when every force is an energy gradient.
Dissipation and pressure work require their own tests and a consistent time scheme.

## 4. Implemented independent geometric oracle

`research/capillary_energy_reference.py` evaluates supplied triangular interfaces.
It does NOT reconstruct an SPH interface and is never called by the production
adapter. For triangle vertices a,b,c, its area is

    A = 0.5 |(b-a) cross (c-a)|,
    dA/da = 0.5 (b-c) cross n,
    dA/db = 0.5 (c-a) cross n,
    dA/dc = 0.5 (a-b) cross n,

where n is the triangle unit normal. The code assembles -gamma*dA/dx on all vertices.
The wetted-area coefficient may be negative, as required by the excess-energy
formula. This is harmless: the omitted solid-gas baseline is constant for the rigid
solid. Degenerate triangles are rejected rather than silently dropping a load.

The tests cover:

- Analytic force versus finite differences of the SAME energy.
- Translation and rotation invariance and zero total force/torque.
- Rigid-body translation and rotation virtual work, including nonzero torque.
- Spherical-cap contact-line uplift converging toward the analytic pressure load
  for 60, 90 and 120 degree caps under geometric refinement.
- Young's angle as an energy minimum along fixed-volume spherical caps.

Marking contact vertices as solid-attached is used ONLY to take a rigid variation
at frozen contact geometry. Permanent attachment would pin the contact line. A real
wetting model needs contact-line motion relative to the solid and must include its
variations. The tests establish force algebra and geometric-reference behavior,
not an SPH contact-angle result, topology handling, or a calibrated capillary model.

## 5. SPH implementation requirements before enabling a replacement

The runtime remains particle-based. A triangle mesh is useful as an independent
oracle but is not being introduced as a persistent interface tracker. The next
SPH operator must construct differentiable liquid-gas and solid-liquid area
estimators from local particle neighborhoods. Necessary requirements are:

1. A common area/energy definition for fluid and solid contributions. It must
   converge to the physical areas; matching only the net force is insufficient.
2. Full derivatives of every position-dependent kernel, normalization, surface
   normal, and solid geometry contribution. Stopping a derivative through a
   reconstructed normal destroys the virtual-work guarantee.
3. A deliberate choice of volume weights. Fixed reference volumes simplify the
   variational force but must be checked under compressibility. If current volumes
   or densities depend on positions, their derivative cannot be omitted while
   claiming a full energy gradient.
4. Euclidean-invariant local geometry. Fixed Cartesian quadrature, clipping,
   hard surface flags or hidden geometry-dependent coefficients need explicit
   examination; conservation cannot simply be assumed for them.
5. Compact smooth support and correct derivative accumulation across neighbor
   changes. Analytic derivatives/automatic differentiation are appropriate for
   verification; production must avoid dense global automatic-differentiation tapes.
6. Derivatives with respect to solid sample positions or rigid pose, not just fluid
   coordinates. At particle overlap, only the physically exposed solid union should
   contribute to wetting; summing buried surfaces double-counts area and adhesion.
7. Replace the existing CSF and prescribed-normal capillary closure when enabling
   the new energy operator. Do not stack independent capillary models. Keep pressure
   and viscous traction separate, and audit their existing reciprocity/time levels.

The concrete SPH area estimator is NOT yet selected or validated. Consequently no
new runtime force is supplied in this update. This is the key remaining modeling
step, not an implementation detail that can safely be guessed from the residual.

## 6. Validation sequence and cost control

Before a fluid trajectory, test the SPH capillary operator on frozen particle
configurations: energy-gradient consistency, translation/rotation invariance,
force/torque reciprocity, planar-interface energy per area and solid-liquid excess
energy per area. Surface tension and contact angle require their own area-energy
convergence tests, not a global force rescaling.

Only after those checks should we run a short drop/plate balance and a fixed
sphere-bridge force/torque comparison against the independent Young-Laplace reference.
Then test moving rigid bodies, split/join events and exposed-surface handling.
Retain the acoustic viscosity correction. Do not request another 20-case sweep or
long equilibration until the local operator checks pass.

## 7. Literature and scope

- Probst & Teschner, *Unified Pressure, Surface Tension and Friction for SPH Fluids*,
  ACM Transactions on Graphics 44(1), article 7; published online 2024.
  https://doi.org/10.1145/3708034
  This paper provides an energy-based particle surface-tension construction that
  includes boundary contributions and discusses momentum conservation. It supports
  the shared-energy principle; our triangle oracle is an independent reference,
  not a reproduction of its complete incompressible solver or friction model.
- Zoeller, Adams & Adami, *A partitioned continuous surface stress model for
  multiphase smoothed particle hydrodynamics*, JCP, 111716 (2023).
  https://doi.org/10.1016/j.jcp.2022.111716
  A relevant alternative for high density ratios and wetting. Its abstract does
  not establish every discrete rigid-body angular-momentum property required here.
- Whebell et al., *Computing sessile droplet shapes on arbitrary surfaces with a
  new pairwise-force smoothed particle hydrodynamics model*, JFM 1019, A37 (2025).
  https://doi.org/10.1017/jfm.2025.10601
  A useful central-force alternative, but contact-angle/surface-tension calibration
  and its focus on equilibrium shapes must not be confused with demonstrated
  conservative moving DEM coupling for this project.

The equations and triangle-force tests above are our derivation/reference work.
They do not establish that any cited paper supplies this project's complete DEM
coupling or that SPHinXsys's pinned curvature operator is variational.
