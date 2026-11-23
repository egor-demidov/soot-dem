> Update: see BULK_PRESSURE_AND_OVERLAP_FIX.md for the subsequent bulk-density derivation and overlap correction. This document records the earlier investigation.

# Fixed-boundary audit and two fresh controls

No force-law change or proposed physical fix is included in this update.

## Equation audit

Reference: Whebell et al., JFM 1019 A37 (2025),
https://doi.org/10.1017/jfm.2025.10601, sections 2.2 and 2.4.

With e_ij = (x_i-x_j)/r, V_i=m_i/rho_i, H=4 dx:

- Eq. 2.4: density rate = rho_i sum_j V_j (v_i-v_j) dot grad_i W.
- Eqs. 2.5-2.6: pair pressure force = -V_i V_j (P_i+P_j) W'(r) e_ij.
- Eq. 2.11: pair viscous force = 10 mu V_i V_j W'(r)/r
  ((v_i-v_j) dot e_ij) e_ij in 3D.
- Eq. 2.7: W' = -210/(pi H^4) q (1-q)^3, q=r/H.
- Eqs. 2.19-2.21: pair capillary force = H s_ij f_ij(q) e_ij;
  polynomial alpha_ff=-23.5 and alpha_fs=-11.
- Section 2.2.4: rigid dummy samples hold rho0 and mass rho0 dx^3;
  their velocities vanish for a fixed sphere and EOS pressure is zero.
- Eqs. 2.15-2.16: midpoint force/density evaluation, endpoint trapezoidal
  position drift. The fixed-body fluid path has this structure.

The current algebra matches these terms. This does NOT establish that every
aspect of sampling, initialization, or numerical implementation is correct,
nor that this finite-resolution curved boundary meets impermeability needs.
Our extra viscous-row/body timestep restrictions are conservative additions.
Moving rigid samples, torque reactions, and overlap masking are extensions;
the original work's fixed-substrate validation does not validate those extensions.

Initialization creates fluid at half-offset global lattice nodes and boundary
samples on a body-centered cubic lattice clipped to each sphere. Both have dx
spacing, but these are independent lattices; there is no common interface gap
or uniform inter-phase sample spacing enforced. This is a discretization choice
that must also be tested, not assumed equivalent to every setup in the paper.
The 14 nm sphere has R/H=2, and the 18 nm droplet has R/H about 2.57. These are
small geometries relative to kernel support. They are useful failure-isolation
tests, not a grid-converged capillary-force benchmark.

## Verified checkpoint result

The supplied aggregate checkpoint at step 40000 has 5163 liquid samples;
1210 are inside the sphere union and 177 penetrate deeper than dx. Minimum
surface gap is -5.94575 nm. Sample 2404 penetrates particle 2 by 2.42657 nm,
with relative outward normal velocity -17.9409 m/s. Force decomposition:

| Component | Outward normal force (pN) |
|---|---:|
| Solid pressure | 401.25264 |
| Solid pair force | -295.56168 |
| Solid viscosity | 52.44534 |
| Liquid pressure | -317.74297 |
| Liquid pair force | 146.22339 |
| Liquid viscosity | 9.50598 |
| Total | -3.87729 |

A separate C++ evaluation of the checkpoint's fluid/rigid state agrees with
Python's summed force and density rate to roundoff. No mechanical spring state
is needed to evaluate these instantaneous fluid terms. The Python audit reads
the checkpoint prefix only and never restores or modifies the mechanical state.
EOS pressure is not total mechanical pressure in this pair-force formulation.

The user's fixed-sphere run reproduces substantial penetration (-5.70 nm at
its maximum-density sample by 0.4 ns) without DEM motion. Investigation therefore
starts with the fixed boundary, not force-lever transfer or moving necks.

## Run two controls

From your simulation directory, after installing this component update:

```bash
ROOT="$HOME/CLionProjects/soot-dem"
python3 "$ROOT/capillary/tools/run_boundary_controls.py" \
  boundary_fixed_sphere.xml \
  --executable "$ROOT/cmake-build-release-gcc-12/pf_soot_dem" \
  --output results-boundary-controls
```

Use your existing boundary_fixed_sphere.xml from the previous test. An identical
reference is bundled at examples/boundary_fixed_sphere.xml. The previous ratio-6
baseline need not be rerun; --include-baseline is optional.

- `no_solid_pf`: sff=0.002, sfs=0. Cohesion and hydrodynamics remain.
- `ratio3`: sff=0.002, sfs=0.006. Same fluid, geometry, sound speed and viscosity.

Removing solid PF removes BOTH adhesion and short-range PF repulsion. It is not
a physical nonwetting material nor an impermeability validation by itself.
Together with the known ratio-6 result, these controls test dependence on the
solid pair strength. They do not uniquely distinguish attraction from repulsion,
nor automatically identify a cure. A pass at ratio 3 would not establish that
the intended contact angle has been recovered.

Runs are fresh (no compressed checkpoint). Each inherits the original end time,
step cap, resolution and output cadence; checkpoint cadence is reduced to 1000.
The script runs the two cases sequentially, reports failed cases, and continues
to the next one. It records exit codes and the last history row in
control_summary.json. It also audits the most recent checkpoint into each
case's boundary_audit.json. Audit time is printed explicitly: if a run fails,
this may precede the final density report/failure snapshot.

Send control_summary.json and both boundary_audit.json files, plus the final
rows of density_extrema.csv. No full trajectory is needed.

## Audit any existing checkpoint without integration

```bash
python3 "$ROOT/capillary/tools/audit_boundary_checkpoint.py" \
  results-sph-levelset/dem_checkpoint.dat --output boundary_audit.json
```

Requires NumPy. The two probe samples are the minimum and maximum density at
that checkpoint; the deepest-penetrating sample can be a different sample.
Normal direction is outward from the nearest sphere. Every geometric liquid
neighbor and unburied solid neighbor inside support is included. The audit
uses stored fluid/boundary data; it does not advance, fix or restart the run.

## Local verification

Three runner/audit tests passed. Each control completed exactly one tiny SI
step (1e-18 s), verifying executable compatibility and reporting only. No
0.4 ns control or long trajectory was run locally. Existing force physics and
executable were reused; the new scripts require no C++ rebuild.
