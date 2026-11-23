# Bulk pressure explanation and overlap quadrature correction

## Meaning of density error

The driver reports max_i abs(rho_i/rho0 - 1) using each liquid sample's evolved
continuity-equation density. rho0 is the input reference liquid density, not a
density recalculated from an initial grid. Laplace pressure and the pair-force
model's bulk cohesive stress can both require positive EOS pressure/compression.
The 30% operating threshold does not distinguish physical from numerical causes.
Treating a positive density offset alone as proof of instability was incorrect.

## Static bulk calculation: no simulation or fitted correction

For an isotropic uniform 3D distribution, the PF virial pressure is

P_pf = (2 pi/3) n^2 integral_0^H r^3 [H sff f_ff(r/H)] dr.

With H=4 dx, m=rho0 dx^3, n=rho/m and alpha_ff=-23.5:

integral_0^1 q^3 f_ff(q) dq = (9+alpha_ff)/168 = -29/336,
P_pf = -C (sff/dx) (rho/rho0)^2,
C = (2 pi/3) 4^5 (29/336) = 185.1046338115.

Balancing the Tait EOS with this negative cohesive stress at zero total bulk
pressure gives (r=rho/rho0):

rho0 c^2/7 (r^7-1) - C (sff/dx) r^2 = 0.

For dx=1.75 nm, rho0=1258 kg/m^3, c=1000 m/s and sff=0.002 N/m:
- predicted bulk ratio r = 1.14215708644;
- EOS pressure = +275.9694 MPa;
- cohesive pressure = -275.9694 MPa;
- total bulk pressure approximately zero.

Observed median ratios in the two supplied controls were 1.14408 and 1.14441.
This close agreement supports cohesive/EOS balance as the dominant bulk offset.
It does NOT predict every local density, curvature stress, boundary effect,
finite-lattice error, or transient. It assumes isotropic homogeneous continuum
sampling and ignores kinetic stress. The finite droplet still has Laplace stress.

Thus the 14% density offset is largely explainable by the chosen PF model and
compressibility, not evidence by itself of a sign bug or density drift. It also
means rho0 is not automatically the equilibrium zero-mechanical-pressure density
of this cohesive liquid. Surface tension, actual equilibrium density, sound speed
and wettability need joint calibration if interpreting the input as a material.
No EOS shift, artificial repulsion, density reset or clipping is added here.

Run the static calculation:

```bash
python3 capillary/tools/pf_bulk_balance.py
```

## Corrected: missing volume quadrature in sphere intersections

The old overlap rule removed every interior site lying inside another sphere.
It deleted the intersection from both body grids. Since the sites sample volume,
a surface-exposure mask is inappropriate.

The new rule assigns each rigid sample a transient weight 1/N(x), where N(x) is
the number of solid spheres containing its point (including its owner). This is
the partition of unity chi_k(x)/sum_j chi_j(x): summed body integrals cover the
solid union once. Each body's interior lattice approximates its own integral.
On coincident grids, duplicated volume samples contribute exactly one copy.
On independently shifted/rotated grids this remains a quadrature approximation.

Weights multiply ALL fluid-solid pair contributions: pressure, viscosity,
capillarity, density rate, pair energy, virial and viscous stability rows.
Equal/opposite forces and common-point torques are retained. Physical liquid
mass and the stored rigid sample mass/density are not changed; weights are
recomputed at force-stage positions. Nonoverlapping bodies follow the existing
weight-1 path. The correction adds no new force profile or fitted coefficient.

This fixes the confirmed volume-deletion defect. It is NOT a validation of
moving overlapping wet particles. Discontinuous cover changes at finite sample
spacing can still introduce force noise, and position-dependent quadrature
weights mean the reported weighted pair energy is not a complete conservative
potential for moving overlaps. Global energy consistency remains unresolved.

## Checkpoints and audits

New solver checkpoints have header PFDEM 3. The audit recognizes version 3 and
uses partition weights; old version 1/2 checkpoints are audited with their old
mask so historical force decompositions retain their meaning. Old checkpoints
with overlapping spheres are refused by the new solver rather than silently
changing their boundary rule on restart. Start fresh for overlapping aggregates.
Old nonoverlapping checkpoints remain readable.

The diagnostic boundary neighbor count now counts all geometric solid neighbors
within support (unweighted); force/density rates include the partition weights.

## Verification / use

Native regressions check overlap coverage (two coincident samples count once,
not zero), pressure/density/viscosity/capillary weighted linearity and reciprocal
force/torque transfer. Existing single-body/nonoverlap and restart checks apply.
Python tests check the independent profile integral and static bulk equation.
No long trajectory is required to reproduce these algebraic/geometric defects.

Reinstall the component and rebuild pf_soot_dem for the overlap correction.
Do not rerun the stationary-sphere controls merely to test this change: their
one-body trajectory should be unchanged. This correction cannot cure the
ratio-6 penetration previously observed in a single sphere.
