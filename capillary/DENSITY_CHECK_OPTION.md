# Optional reference-density stopping check

Default behavior is unchanged: stop when max |rho/rho0-1| > 0.3 at a logging
check. To disable ONLY that screen, either add density_check="false" to the
<sph> section of the original <simulation> XML, or to <run> in modern <pf_dem>
XML. Alternatively pass --disable-density-check to tools/run_pf_dem.py or to
the C++ executable. Python translates the XML runtime setting into that flag;
the numeric prepared-input format and model signature are unchanged.

Density ratios/pressures are still calculated and reported. The Mach 0.2 screen,
nonpositive/nonfinite state validation, EOS overflow and timestep checks remain.
This switch does not assert validity of a state or repair penetration. It does
not change forces, viscosity or density integration.

## Continue the original levelset trajectory

After installing this component update and rebuilding your actual executable:

```bash
ROOT="$HOME/CLionProjects/soot-dem"
cmake --build "$ROOT/cmake-build-release-gcc-12" --target pf_soot_dem -j 6
python3 "$ROOT/capillary/tools/run_pf_dem.py" capillary_levelset_small.xml \
  --executable "$ROOT/cmake-build-release-gcc-12/pf_soot_dem" \
  --resume "$PWD/results-sph-levelset/dem_checkpoint.dat" \
  --disable-density-check --legacy-overlap-mask \
  --end-time 1e-9 --output results-sph-levelset-continued
```

The sample target is 1 ns; change --end-time for a different physical time.
The source checkpoint and output directory are not overwritten. Restart begins
at the last saved checkpoint, not necessarily the failing log step.

## Preserve the trajectory's original overlap rule

--legacy-overlap-mask explicitly allows a v1/v2 checkpoint with overlapping
spheres to continue with the old mask. This compatibility mode reproduces the
old rule, including its known volume-deletion defect; it does not silently
switch an old trajectory to the corrected partition. It is available only with
--resume, and has no effect on a checkpoint already using the corrected rule.

New checkpoints use PFDEM 4 and store which overlap rule is active. A later
restart preserves it without requiring the legacy flag again. The audit tool
understands both stored rules. Unflagged old overlapping checkpoints continue
to be rejected. Existing PFDEM 3 partition checkpoints remain supported.

## Verification

Tests check that disabling the XML screen leaves the model signature and numeric
prepared input unchanged, CLI flags propagate, the native solver stops at a
synthetic 31% density deviation by default, continues for one tiny step when
disabled, and still stops for excessive Mach number. Static legacy-mask tests
verify explicit old-checkpoint continuation and persisted PFDEM 4 rule state.
No long trajectory is run locally.

The uploaded real checkpoint was not advanced locally: its soot-dem mechanical
source fingerprint differs from this local checkout. That guard remains intact.
Use your original local soot-dem source revision when rebuilding to resume it.
