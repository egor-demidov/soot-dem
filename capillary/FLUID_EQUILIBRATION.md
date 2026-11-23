# Initial fluid equilibration

Add an optional duration in seconds to the existing soot-dem XML:

```xml
<let id="fluid_equilibration_time" type="real">1e-9</let>
```

Alternatively put `fluid_equilibration_time="1e-9"` on the existing `<sph>`
section. Use one location only. For the modern `<pf_dem>` input, put the
attribute on `<run>` instead. The default is zero, preserving immediate coupled
motion. Both spherical-droplet and levelset initialization support this option.

During `[0, fluid_equilibration_time)`, all DEM particle centers and orientations
are fixed. Their translation and spin are zero at the fluid boundary, while
fluid positions, velocities and density evolve with the usual operators. Native
soot-dem mechanical force evaluations and history commits are paused: necks are
retained, and tangential/rolling/torsion spring histories do not accumulate.
Fluid reactions on the clamped particles are assigned to the support impulse
ledger. No fluid load is accumulated for an instantaneous kick at release.

At the release time, particles recover their original fixed/free flags and their
configured initial linear/angular velocities. Permanently fixed particles stay
fixed. Initially stationary particles start coupled motion from rest under the
current fluid and mechanical loads. A timestep is shortened to hit this release
time exactly, and rates are recomputed after release.

The equilibration period **is included** in the total simulation clock and XML
`t_tot`/`end_time`. For example, `t_tot=5e-9` with equilibration `1e-9` gives 1 ns
of fluid-only evolution followed by 4 ns of coupled evolution. A run ending
before release is a fluid-only run and can be continued later. Existing step
cadences remain active during equilibration; if output is enabled, release also
triggers a frame and history record. History adds `fluid_equilibrating` (0/1)
and `coupled_time` (seconds since the most recent release, zero during a hold).

The XML hold starts on a **fresh trajectory**. Checkpoints preserve its
progress, release time, original particle motion and fixity. The XML runner
uses the checkpoint clock on resume. A restart during equilibration completes
only the remaining hold. A restart after release continues the coupled run.
For a new hold after adding liquid, see [PF_CONDENSATION_RESTART.md](PF_CONDENSATION_RESTART.md).
Checkpoints from before the equilibration feature start in the coupled phase.

The Python XML runner forwards the duration as `--fluid-equilibration-time`
to the native executable. This is a runtime control, so prepared-input geometry
signatures are unchanged. New checkpoints use PFDEM version 5; use this source
version to resume them. Readers retain versions 1–4 with their original overlap
compatibility rules.

This provides a fixed-duration relaxation stage, not an automatic equilibrium
criterion. Inspect fluid velocities and loads to decide whether a chosen hold
is long enough. It does not relax the dry DEM contact/neck geometry.

## Short verification

`pf_dem_test` checks frozen poses, fluid evolution, release of initial velocities
and spin, permanently fixed bodies, support momentum accounting, exact mid-hold
restart and continued DEM motion. `pf_soot_dem_test` additionally checks paused
native spring histories, preserved necks, native checkpoint restoration and
resumption of history updates after release. Python tests cover both XML layouts,
default/invalid durations, argument forwarding, native release timing and restart
conflict handling. No long relaxation trajectories are required by these tests.
