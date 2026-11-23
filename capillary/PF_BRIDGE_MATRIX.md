# Force calibration with a 2×2 bridge matrix

The matrix uses fixed equal spheres, radii of 21 and 31.5 nm, and
surface-to-surface gaps of 7 and 10.5 nm. Every radius is paired with both
gaps. All four runs use the same fluid parameters, including `dx=1.75 nm`
and `sff=0.002 N/m`. The default `sfs=0.008 N/m` characterizes the current
bridge runs. The seed droplet radius is 5/6 of the sphere radius and the
Young–Laplace calculation uses the measured late liquid volume of each
case.

## Run from the simulations directory

```bash
TOOLS="$HOME/CLionProjects/soot-dem/capillary/tools"
cd /media/egor/Data/sph-dem-simulations

python3 "$TOOLS/pf_bridge_matrix.py" prepare bridge-matrix-sfs008 --sfs 0.008
python3 "$TOOLS/pf_bridge_matrix.py" run bridge-matrix-sfs008 --executable pf_soot_dem
python3 "$TOOLS/pf_bridge_matrix.py" analyze bridge-matrix-sfs008 \
    --gamma 0.08781 > bridge-matrix-sfs008/force_angle.json
```

`prepare` creates a fresh directory with four editable XML files.
`run` skips completed cases. To reuse a previous run, link its directory
to the matching case name under the matrix directory before running
(e.g. `r31.5-g10.5`). Analysis checks the sphere radius, gap, fixity, and
fluid parameters across the matrix. Inspect the final VTPs for bridge
attachment and liquid penetration, especially near rupture or strong wetting.

Analysis averages the axial force and liquid volume over the final
quarter of recorded time. It requires completed runs with at least five
late records and screens late force variation, force imbalance, volume
variation, and Mach number. It uses signed attraction, so repulsive
bridges are analyzed as negative forces.

Each case gets a force-equivalent contact angle when the reference curve
crosses its SPH force. The common angle minimizes the mean squared force
error normalized by `gamma × sphere radius`. A separate fit holds out
`r31.5-g7` and predicts its force from the other three cases. The angle
is an effective calibration for this discretization and force model.
The default search is 15–165 degrees. Use `--theta-range LOW HIGH` to narrow
or expand the scan within (0, 180).

## Obtuse analytical angles

`young_laplace.bridge_reference(volume, gap, theta_deg, radius, sigma)`
now accepts angles strictly between 0 and 180 degrees through the liquid,
including 90 degrees. Its arclength formulation permits axial overhangs,
and it retains the sign of the force. Negative attraction means
repulsion. It searches for a simple symmetric meniscus contacting the
facing hemispheres of equal spheres at a positive gap. Some volumes and
gaps have no solution on this branch. A converged solution does not
establish stability.

The analytical angle is separate from the SPH XML. To characterize a
different fluid–solid interaction, prepare another matrix with a different
`--sfs`. Use the appropriate independently calibrated surface tension.
