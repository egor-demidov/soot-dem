> Update: see BULK_PRESSURE_AND_OVERLAP_FIX.md for the subsequent bulk-density derivation and overlap correction. This document records the earlier investigation.

# Targeted bug review after fixed-sphere controls

## Fixed: empty boundary neighborhood in the checkpoint audit

A minimum/maximum-density liquid sample can have no solid samples within H.
The old audit called min() on that empty distance array and crashed, losing the
whole report after the solver had completed successfully. It now writes JSON
null for nearest_boundary_sample_distance_dx, zero solid counts and zero solid
force vectors. Empty support is not evidence of an invalid fluid state.

Regression coverage also handles a liquid sample exactly at the center of a
sphere: radial normal components are undefined and represented by null, rather
than a divide-by-zero producing NaN. The audit now rejects distinct coincident
active samples, consistent with the solver, instead of silently omitting them.
All these changes are diagnostic only; no trajectory or checkpoint changes.

Three Python tests pass, including empty-neighbor and undefined-normal cases.
No new trajectory benchmark was run for this update.

## Confirmed separate solver flaw: volume mask in overlapping spheres

InitialGeometry fills sphere interiors with a volume lattice. System::stage
then suppresses every site's interactions if its position is inside another
sphere. In an intersection region, sites belonging to both spheres are
suppressed. This creates missing volume quadrature, not a single representation
of the solid union. The same rule was used in the audit, so it faithfully
reproduced that solver behavior rather than detecting the geometry flaw.

A static C++ fixture confirms it: one liquid sample outside two overlapping
spheres, two interior sites in their intersection, both in liquid kernel
support. There are two unmasked edges and zero active edges after stage masking.
This proves the all-sides deletion. It does not quantify the effect on a real
trajectory. It is absent in a one-sphere control and cannot explain that case.

The production overlap treatment has NOT been changed in this diagnostic fix.
It needs a volume partition/ownership discretization with force/torque transfer,
not a surface-exposure test applied to interior quadrature. Choosing ownership
also affects which DEM body receives loads in the overlap; that choice and
switching behavior need to be explicit. This is separate from the current
stationary single-sphere investigation.

## What the submitted controls establish

Both simulations reached 0.4 ns successfully. no_solid_pf's audit alone failed.
ratio3 reports 0 of 4124 liquid samples inside the sphere; minimum surface gap
0.66838 nm. Maximum density error remains 15.927%, so this is not a calibrated
weak-compressibility/contact-angle validation. no_solid_pf's maximum density
error is 22.034% and max speed 37.70 m/s; its penetration counts require the
repaired audit. The ratio-6 reference showed substantial penetration.

These results refute an inference that the stationary boundary always fails in
the same way. They demonstrate strong dependence on solid pair strength, but
do not by themselves identify a single implementation error causing ratio-6
penetration. The force signs, volume factors, kernel derivative, EOS, viscosity
and continuity algebra checked so far match the reference paper. The initial
sampling uses two independent lattices, and the small radii relative to H remain
unvalidated discretization choices, not confirmed algebraic bugs.

## Recover the missing report without rerunning simulations

After installing this scripts-only update:

```bash
ROOT="$HOME/CLionProjects/soot-dem"
python3 "$ROOT/capillary/tools/audit_boundary_checkpoint.py" \
  results-boundary-controls/no_solid_pf/dem_checkpoint.dat \
  --output results-boundary-controls/no_solid_pf/boundary_audit.json
```

The old control_summary.json retains the old audit_exit=1; this standalone audit
command intentionally does not rewrite historical runner exit status. Upload the
new boundary_audit.json. No C++ rebuild or solver run is needed for this recovery.
