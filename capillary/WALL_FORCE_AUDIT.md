# Continued fixed-boundary investigation: static wall-force balance

This is a mathematical diagnosis, not a new wall force, density correction,
contact-angle calibration or trajectory benchmark.

## Checked against authors' code

The paper's HTML link omits an underscore and failed. The accessible public
reference is https://github.com/rwhebell/Whebell2024_SessileDroplets,
inspected at commit 1164577b1b879fb9f5d3e971caa5bf89c89afa0b.

`deltaplusSPH.jl` functions stateEquation!, calcInteractions! and the midpoint
update were inspected directly, alongside surfaceTension.jl and kernels.jl.
The checked EOS, rho_i V_j continuity rate, (P_i+P_j) pressure term, 10 mu radial
viscosity and H s/m_i pair acceleration agree with our fixed-boundary equations.
Our undirected pair implementation assembles the equal/opposite forces once.
The reference initializes fixed-ghost pressure to zero and updates EOS pressure
only for liquid indices; our rho0 ghosts likewise have zero Tait pressure.
No missing pressure normalization or pair-strength factor was found in those
terms. This is narrower than claiming complete equivalence of implementations.

The authors' continuum model and calibration concern fixed substrates. Our
moving-body reaction, overlap partition and rigid histories are extensions.
No Julia simulations, or other long trajectories, were run for this audit.

## Force at an ideal plane

Assume the solid occupies z<0 and the liquid fills z>0. Both distributions are
homogeneous in their respective half spaces. Dummy number density is dx^-3;
liquid number density is r dx^-3, r=rho/rho0. Pressure is P in liquid, zero in
solid, and velocities vanish. Evaluate normal force on one liquid sample at
z=0, positive outward from the solid.

For the published polynomials:
  integral_0^1 q^2 f_alpha(q) dq = (10+alpha)/105.

Thus I_ff=-13.5/105 and I_fs=-1/105. The combined pressure term of both half
spaces is NOT the solid pressure term alone: their sum is
  F_pressure = -(3/8) P dx^2/r.

Pair-force contributions are
  F_ff = pi*4^4*(13.5/105)*sff*dx*r,
  F_fs = -pi*4^4/105*sfs*dx.

Even though an individual fluid-solid pair repels at short distance, the
half-space integral of its pair-force profile at the plane is attractive.
Consequently, increasing sfs increases attraction and repulsion together but
does not guarantee stronger NET exclusion at the geometric wall.

At the bulk-balance density derived in BULK_PRESSURE_AND_OVERLAP_FIX.md,
P=185.1046338115*(sff/dx)*r^2, giving
  F_normal = (pi*4^4/105)*sff*dx*((71/16)*r - sfs/sff).

For the supplied parameters r=1.14215708644:

| Contribution on one liquid sample | Ratio 3 | Ratio 6 |
|---|---:|---:|
| Combined pressure, pN | -277.48687 | -277.48687 |
| Fluid-fluid pair force, pN | +413.35975 | +413.35975 |
| Fluid-solid pair force, pN | -80.42477 | -160.84954 |
| Total, pN | +55.44811 | -24.97666 |

The zero at the geometric plane is sfs/sff=(71/16)*r=5.06832207108.
This is a plane/uniform-density estimate, NOT a universal allowed-strength
limit, full potential barrier, discrete curved-wall stability criterion or
prediction of the eventual penetration depth. Local density changes, curvature,
fluid distribution, velocity, finite quadrature and transients alter it.
Nevertheless its change of sign is consistent with the observed ratio-3 and
ratio-6 results. Strong adhesion can defeat effective exclusion in this dummy
boundary model without a force-sign programming error.

## Practical conclusion

Two distinct findings must remain separate:
1. The overlapping-volume deletion was a confirmed implementation defect and
   was corrected in the preceding update.
2. The ratio-6 single-sphere result is consistent with a limitation of the
   chosen finite-resolution dummy-particle wetting model. Algebra checked against
   the authors' implementation has not exposed an additional sign/factor bug.

Use ratio 3 only as the empirically nonpenetrating diagnostic setting observed
so far. It does NOT establish the intended contact angle. Supporting stronger
wetting with reliable impermeability requires a separately justified boundary
formulation or resolving/calibrating this one; adding an arbitrary repulsive
term would change the model and its force/contact-angle calibration.

The driver message now explicitly calls the 30% check a reference-density
operating screen. The historical max_density_error CSV column is retained for
compatibility, but denotes deviation from rho0 rather than numerical error.
No threshold or force is changed in this update.

Run the STATIC calculations (no particle dynamics):

```bash
python3 capillary/tools/pf_bulk_balance.py
python3 capillary/tools/pf_wall_balance.py --sfs 0.006
python3 capillary/tools/pf_wall_balance.py --sfs 0.012
```

Tests independently integrate the polynomial profiles and kernel gradient,
verify the bulk-balanced threshold, and check the opposing force signs.
