# Reviewed dependency

SPHinXsys v1.0-beta.08, commit 319520e02f631847f85160c0e8b55648c395eafa.
This version contains the explicit contact-angle class documented in example 18.
The current master inspected during setup (1612330e44c8947406f24b6756449e80dcd7158e)
contains a newer surface-stress implementation and no SurfaceNormWithWall class.
We deliberately do not combine methods from these two versions.

Reviewed files in the pinned dependency:
- tests/2d_examples/test_2d_wetting_effects/src/wetting.cpp
- src/shared/particle_dynamics/fluid_dynamics/fluid_surface_complex.h
- src/shared/particle_dynamics/fluid_dynamics/fluid_surface_inner.hpp
- src/shared/particle_dynamics/solid_dynamics/fluid_structure_interaction.h

SurfaceNormWithWall rotates a liquid interface normal toward the input angle.
SurfaceTensionAccelerationInner adds a CSF-like curvature force to liquid acc_prior.
The FSI class computes pressure forces using pressure, a Riemann dissipation term,
and acc_prior in wall pressure extrapolation. Capillary acceleration can therefore
enter that extrapolation indirectly. There is no separate solid wetting-reaction
class in these reviewed routines. This audit does NOT prove that capillary force
is absent or correct. The bridge-force and supported momentum tests assess it.

Adapter additions are explicit case geometry, diagnostic output, a capillary time
step restriction, and experimental dry-sphere integration. No global conservation
projection, contact-line force estimate, mass restoration, or upstream physics
patch is used. The moving case uses a first-order split symplectic Euler update
with sphere inertia and quaternion tracking; it is an audit, not validated DEM
coupling. Impulse sampling and splitting errors require time refinement too.
