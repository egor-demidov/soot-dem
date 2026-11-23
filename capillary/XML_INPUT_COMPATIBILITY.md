# Original soot-dem XML layout for PF-SPH coupling

`tools/run_pf_dem.py` accepts both the existing `<pf_dem>` files and the original
`<simulation><let>...<include>...<initialization>...` layout. The latter selects
the actual soot-dem force backend; run it with **pf_soot_dem**, not pf_dem.

## Updated original inputs

Copy **all four** files from `inputs/` beside your existing `aggregate_bank/`:

```bash
cp capillary/inputs/capillary*.xml .
python3 capillary/tools/run_pf_dem.py capillary_levelset.xml \
  --executable "$PWD/build/pf_soot_dem" --output results-sph-levelset
# Alternative coating, using the original droplet coordinates:
python3 capillary/tools/run_pf_dem.py capillary_dem.xml \
  --executable "$PWD/build/pf_soot_dem" --output results-sph-droplet
```

After extracting this component update, use its `tools/install_into_soot_dem.py`
installer to update `capillary/` without replacing your soot-dem root CMake.
There is no C++ model change in this input-format update. Existing executables
from the preceding soot-dem integration can be reused.

Paths are relative to the XML file containing the path. Nested include paths
are relative to their own include file. The example aggregate remains
`aggregate_bank/aggregate_1/final_aggregate.vtk`; that user file is not bundled.
Original soot-dem VTK positions are in primary-radius units, now explicitly
selected by `length_unit="primary_radius"`. SI or nm coordinate files should
instead use `length_unit="m"` or `"nm"`. Imported velocities are ignored.

| Original field | Current meaning |
| --- | --- |
| `r_part`, `rho` | Primary radius in m and solid density in kg/m³ |
| `liquid_density`, `liquid_dynamic_viscosity` | Fluid density in kg/m³ and dynamic viscosity in Pa·s |
| `dt` | Maximum DEM integration timestep in seconds; adaptive stability limits may reduce it |
| `t_tot` | End time in seconds |
| `initialization/aggregate` | Existing aggregate geometry; automatic touching-particle neck graph by default |
| `initialization/liquid mode="droplets"` | Union of spheres with solid primary-particle interiors excluded |
| `initialization/liquid mode="levelset"` | Sum of erf profiles, thresholded outside primary-particle interiors |
| `liquid rm` | Profile radius in the specified liquid length unit, independent of solid radius; omitted/zero uses each primary radius |
| `sph particles_per_primary_diameter` | SPH sample resolution, replacing LB lattice resolution |
| `domain_margin` | Per-axis geometric sampling margin; no periodic images or wall are created by this box |
| `output ..._every_sph_steps` | Accepted coupled integration-step intervals, not fixed LB timesteps |

The levelset copy preserves the latest uploaded settings: alpha 0.5, threshold
0.5, 14 nm radius, 60° target, 10 nm margins. The droplet copy preserves its
original center (-140.5052, 90.9034, 143.8088) nm, radius 60 nm, 5° target,
and 30 nm margins. Both preserve 14 nm primary radius, density 1700 kg/m³,
fluid density 1258 kg/m³, viscosity 0.001 Pa·s, dt 5e-13 s and end time 5e-8 s.

## Calibration and neck parameters

`liquid_surface_tension` and `liquid_contact_angle_deg` are retained as
**physical targets**, written to `resolved_input.json`. They are not silently
converted into an asserted calibration. This pair-force formulation needs
explicit `sph sff` and `sph sfs`. The supplied sff = 0.063/30.96 and sfs/sff = 3
are provisional starting estimates; neither example is a calibrated 5° or
60° case. `sound_speed` and particle spacing are SI quantities. Change the
explicit force strengths after calibration, retaining the target metadata.

`capillary_contacts.xml` preserves the recovered original dry-contact and
van der Waals coefficients. That old include contained **no neck stiffness
coefficients**. `capillary_necks.xml` therefore supplies clearly marked NEW
provisional coefficients, matching the corresponding contact stiffness and
damping, with a 1 nm automatic neck tolerance. Review these for your material.
Zero yield thresholds disable neck breakage. Necks contribute only DEM forces,
never fluid geometry or volume. Use `<necks mode="none"/>` to disable them, or
`<necks mode="file" file="necks.vtk"/>` for an explicit saved topology.

## Output and legacy compatibility

The driver currently shares one fluid/particle dump cadence and one DEM/log
record cadence. Set the corresponding pairs equal. Zero fluid/particle cadence
disables point-cloud dumps; log cadence must remain positive. Checkpoint cadence
is independent and may be zero. Surface meshes are produced **offline** by the
existing surface exporter; the surface setting records the intended cadence and
must match fluid dumps or be zero. It does not enable automatic online meshing.

Legacy `lb_cells_per_particle_diameter`, `..._every_lb_steps`,
`dem_every_dem_steps` and the VOF `<lattice>` section are recognized for easier
migration. Cell count becomes SPH sample count; the old subcell VOF quadrature
is unused. `r_verlet` and `neighbor_update_period` are recognized as obsolete,
with explanations in resolved metadata. Unknown fields fail rather than silently
changing a simulation. Prefer the new `<sph>` and output spellings in these copies.

## Minimal verification for this update

Six XML compatibility tests and twelve existing DEM/runner/export tests passed.
Both coating modes completed a single native SI step using three synthetic
particles and reduced sample counts (16 and 80 liquid samples). A coarser
levelset grid was correctly rejected as unresolved/empty before using the finer
grid. No long trajectory, original aggregate run, or physical calibration was
performed as part of this input change.
