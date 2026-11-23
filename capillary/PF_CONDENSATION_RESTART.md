# Repeated condensate additions during soot restructuring

Start with the usual levelset initialization and an initial fluid hold.
For the legacy XML use a top-level let parameter.

```xml
<let id="fluid_equilibration_time" type="real">2e-10</let>
```

For modern `<pf_dem>` XML use `fluid_equilibration_time="2e-10"` on
`<run>`. With either XML dialect, this example holds for 0.2 ns and then
lets the aggregate move until the first checkpoint at 1 ns.

From the simulations directory, adjust the XML and executable paths as needed.

```bash
TOOLS="$HOME/CLionProjects/soot-dem/capillary/tools"
XML=capillary_levelset.xml
EXE=pf_soot_dem

python3 "$TOOLS/run_pf_dem.py" "$XML" --executable "$EXE" \
    --end-time 1e-9 --output dose-00

python3 "$TOOLS/run_pf_dem.py" "$XML" --executable "$EXE" \
    --resume dose-00/dem_checkpoint.dat \
    --add-liquid-volume 1e-24 --post-addition-equilibration-time 2e-10 \
    --end-time 2e-9 --output dose-01

python3 "$TOOLS/run_pf_dem.py" "$XML" --executable "$EXE" \
    --resume dose-01/dem_checkpoint.dat \
    --add-liquid-volume 1e-24 --post-addition-equilibration-time 2e-10 \
    --end-time 3e-9 --output dose-02
```

The dose is in m³, so `1e-24` is 1000 nm³. Use a fresh output directory
and the original XML for each stage. `--end-time` is absolute, and the
post-addition hold is a duration from the checkpoint time. Set the end
time beyond release to get a coupled relaxation interval in the same run.

A dose is a one-time restart operation. To continue an interrupted stage,
use its latest checkpoint without either addition flag. A checkpoint
saved during a hold resumes that hold and releases at its stored time.
A new dose requires the preceding hold to have ended. Liquid–solid
coupling runs during the hold. The aggregate poses are fixed and the
mechanical contact and neck histories are paused. The previous fixity
flags, translational velocities and angular momenta are restored at release.

## Stop after kinetic energy decays

Set a fraction strictly between 0 and 1. This example stops at 1% of
the running peak and requires the ratio to stay at or below that target
for 0.1 ns.

```xml
<let id="ke_target_fraction" type="real">0.01</let>
<let id="ke_hold_time" type="real">1e-10</let>
```

For modern XML, put `ke_target_fraction="0.01" ke_hold_time="1e-10"`
on `<run>`. The runner also accepts `--ke-target-fraction 0.01
--ke-hold-time 1e-10`, which override the corresponding XML values.
Without a hold time, it stops at the first threshold crossing.

The reference is the highest total kinetic energy observed in the current
invocation, counting only the coupled phase. The energy is the lab-frame
sum of liquid translation and solid translation and rotation, as in
`dem_history.csv`. The criterion is evaluated each integration step.
A restart starts a new peak and hold timer. Sustained bulk drift can keep
the ratio above the target.

With the option enabled, the legacy `t_tot` (modern `end_time`) is
the absolute time cap. A command-line `--end-time` can shorten it.
An early stop writes the final history, checkpoint, and enabled VTP
outputs and reports `termination_reason: kinetic_energy` in
`dem_run.json`. Reaching the cap reports `time_limit`. The existing
`completed_time` field indicates whether the time cap was reached.

## Volume and surface rule

Select the placement rule for each dose with `--liquid-addition-mode`.

| Mode | Placement |
| --- | --- |
| `liquid-vapor` (default) | Fills sites closest to the liquid-only Wendland volume-fraction contour at 0.5 first. |
| `uniform-liquid-vapor` | Spreads sites across the same eligible liquid–vapor interface. |
| `uniform-solid-vapor` | Spreads sites across exposed portions of the solid spheres. |
| `curvature-based` | Samples liquid–vapor and exposed solid–vapor sites together with Kelvin weights. |

The liquid–vapor modes use lattice sites 0.85–1.45 `dx` from the nearest
saved liquid particle, at least 0.5 `dx` from solid, and on the vapor side
of the contour. The solid–vapor mode samples the shell 0.5–1.5 `dx` from
the nearest solid surface, with no saved liquid particle within
1.45 `dx` and liquid kernel occupancy at most 0.1. These local criteria
also allow gas-filled pockets. Both uniform modes use deterministic
spatial stratification over eligible lattice sites. Small doses produce
discrete patches. The uniform modes use one layer of the saved geometry.

Both `liquid-vapor` and `curvature-based` can add multiple layers in one
dose. When a layer cannot fit the remaining volume, the code fills its
eligible sites at nominal mass and rebuilds the surface using the newly
placed particles. All layers are placed at the checkpoint time, before
the post-addition hold or any SPH evolution. If a later layer runs out
of eligible sites or the sample limit is reached before the target
volume, the entire dose fails.

## Curvature-based dosing

Add these flags to a dose command to use Kelvin weighting.

```bash
--liquid-addition-mode curvature-based --saturation-ratio 1.2 \
    --kelvin-length 1e-9 --liquid-addition-seed 42
```

The saturation ratio is dimensionless and positive. The Kelvin length
is in metres and can be zero to give equal weights when saturation is
above one. At each site, the code uses the positive part of
`saturation_ratio - exp(0.5 * kelvin_length * mean_curvature)`. The mean
curvature is the average of the principal curvatures, positive on a
convex surface viewed from the vapor. Liquid–vapor curvature comes from
derivatives at the nearby 0.5 liquid volume-fraction contour. Exposed
solid–vapor curvature is `1/R` for the nearest exposed sphere of radius
`R`. The liquid interface must resolve at the chosen `dx` for its sites
to participate.

Each successive draw is proportional to the remaining sites' weights.
Sampling is without replacement on the lattice. The seed defaults to
0 and can be changed between doses. The solver recomputes curvatures and
weights for each layer. If no eligible site has a positive weight, the
entire dose fails. The report gives totals across layers of resolved
sites, positive-weight sites, and inserted particles for each interface.

## Volume accounting

Liquid–vapor additions inherit locally interpolated liquid density and
velocity. Solid–vapor additions start at the reference liquid density
and the local rigid-body velocity of the nearest solid. New masses are
adjusted together in the final layer to meet the prescribed initial
`sum(m_i/rho_i)` increment. This adjustment is limited to 20% of nominal
particle mass and reported as `added_particle_mass_scale`. If the
remaining volume is too small to resolve, or a single-layer mode
exceeds its available sites, the run exits before writing a new
checkpoint.

`liquid_addition.json` records the placement mode, layer count, eligible site count,
requested and realized insertion volume, added mass, particle count
and release time. The added mass is
then conserved during the stage. The geometric liquid volume can change
as the densities relax. New VTPs include a `Mass` field, which the
updated reconstruction reader uses for the slightly different masses
of the inserted particles. The volume step is a prescribed addition of
condensate. It does not model vapour transport or condensation kinetics.
