# Surface tension and contact angle with XML and VTP

Use the rebuilt `pf_soot_dem` or `pf_dem` with `tools/run_pf_dem.py`.
A calibration XML contains one `<droplet>` and no DEM particle declarations.
Adding `<wall z="0" half_width="42e-9"/>` cuts the initial droplet at the wall
and creates a static plate. The plate is centered laterally under the initial
droplet and sampled down to the full kernel support. The wall is represented
by the same fixed-density boundary samples used for the spherical DEM
boundaries. The sessile run stops if liquid penetrates more than half a
sample spacing into the wall or comes within one kernel support of the plate
edge near the wall. Increase `half_width` if the liquid spreads to the edge.

The runner writes its usual `liquid_*.vtp`, `dem_history.csv` and
`dem_checkpoint.dat`. The sessile case also writes `wall.vtp` at the
physical wall plane. Use `--resume` with the same XML geometry and
`--end-time` to extend relaxation. Use a fresh output directory.

## Build and run

From the soot-dem root, after installing the updated component, rebuild
`pf_soot_dem` in your existing build directory. These examples use
`cmake-build-release-gcc-12`.

```bash
cmake --build cmake-build-release-gcc-12 --target pf_soot_dem -j 6
C=capillary
EXE=cmake-build-release-gcc-12/pf_soot_dem
```

The XML files use SI units and example fluid parameters. Match `dx`, `rho0`,
`sound_speed`, `viscosity` and `sff` to the aggregate run. The example end
time is 2 ns. Inspect the late-time diagnostics to decide how long to relax.

## Surface tension

```bash
for r in 21nm 24p5nm 28nm; do
    python3 "$C/tools/run_pf_dem.py" "$C/inputs/calibration_gamma_$r.xml" \
        --executable "$EXE" --output "cal-gamma-$r" || break
done

python3 "$C/tools/pf_calibrate_xml.py" gamma \
    cal-gamma-21nm cal-gamma-24p5nm cal-gamma-28nm > gamma.json
```

The free droplets report volume-equivalent radius and total conservative
core virial pressure in `dem_history.csv`. The pressure includes the EOS and
cohesive pair forces. The analysis averages the final quarter of each history
and fits `P_core = offset + 2 gamma/R` using the relaxed radii. In SI, pressure
is in Pa and the slope is in N/m. The intercept, residual, Mach number and
pressure drift help assess the fit. The JSON lists failed screening checks.
If the core is too small, increase the droplet radii at the same `dx`. A
typical starting point is a smallest radius of at least `12 dx`.
Repeat the radius sweep after changing `sff`.

## Contact angle

After choosing `sff`, vary `sfs` in copies of the sessile XML and run each
with its own output directory. The input starts as a hemisphere at the
wall. `sfs` controls wetting, and the angle is measured after relaxation.

```bash
python3 "$C/tools/run_pf_dem.py" "$C/inputs/calibration_sessile.xml" \
    --executable "$EXE" --output cal-angle
python3 "$C/tools/pf_calibrate_xml.py" angle cal-angle \
    --mesh cal-angle/interface.vtp > angle.json
```

The analysis uses the packages in `requirements.txt` (NumPy, SciPy and
scikit-image). It reconstructs a liquid density isosurface from the VTP
points and fits a volume-constrained spherical cap at the known wall plane.
It reports the angle through the liquid, fit error, sensitivity to excluding
points near the wall, and changes over the last three snapshots. Load
`interface.vtp` with `wall.vtp` in ParaView for a visual check. The
reconstructed mesh is an offline output for the last frame. The simulation
VTP series contains the raw samples.

A detached or unresolved drop can yield an apparent angle that fails the
fit checks. Use an attached, relaxed drop and compare results at different
spacings if transferring parameters to another resolution. A flat-wall
angle calibrates this boundary model. The angle and capillary force on
curved DEM particles still depend on how well their geometry is resolved.
