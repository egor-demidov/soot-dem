#!/usr/bin/env python3
"""Run two gas-free PF-SPH trials. This is calibration evidence, not a pass/fail physics certificate."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import zipfile


def analyze(folder):
    meta = json.loads((folder/'pf_run.json').read_text())
    with (folder/'pf_history.csv').open() as f:
        rows = list(csv.DictReader(f))
    first, last = rows[0], rows[-1]
    tail = [r for r in rows if float(r['time']) >= .8*meta['time']]
    def finite_values(name):
        return [float(r[name]) for r in tail if math.isfinite(float(r[name]))]
    pressures = finite_values('virial_pressure')
    result = dict(metadata=meta, physics_status='NOT_VALIDATED',
                  mass_relative_change=float(last['mass'])/float(first['mass'])-1,
                  volume_relative_change=float(last['volume'])/float(first['volume'])-1,
                  max_mach=max(float(r['max_speed']) for r in rows)/meta['sound_speed'],
                  max_density_error=max(float(r['max_density_error']) for r in rows),
                  max_linear_error=max(float(r['linear_error']) for r in rows),
                  max_angular_error=max(float(r['angular_error']) for r in rows),
                  tail_max_speed=max(float(r['max_speed']) for r in tail),
                  equilibrium_established=False)
    if pressures:
        result['tail_virial_pressure_mean'] = sum(pressures)/len(pressures)
        result['tail_virial_pressure_range'] = max(pressures)-min(pressures)
        # One-radius estimate only: use multiple radii for a Laplace validation later.
        estimates = [float(r['virial_pressure'])*float(r['volume_radius'])/2 for r in tail
                     if math.isfinite(float(r['virial_pressure']))]
        result['apparent_surface_tension_one_radius'] = sum(estimates)/len(estimates)
        result['paper_surface_tension_mapping_comparison_only'] = 30.96*meta['sff']
    try:
        from analyze import exterior_sphere_fit
        with (folder/'pf_final_liquid.csv').open() as f:
            result['final_shape'] = exterior_sphere_fit(list(csv.DictReader(f)), meta['dx'],
                                                       'sessile' if meta['case']=='sessile' else 'laplace')
    except ImportError:
        result['final_shape'] = {'status': 'install numpy and scipy for optional shape fit'}
    (folder/'pf_analysis.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', default='pf_capillary_validation')
    p.add_argument('--output', type=Path, default=Path('results-pf-first'))
    p.add_argument('--end-time', type=float, default=5.)
    p.add_argument('--threads', type=int, default=6)
    p.add_argument('--sff', type=float, default=.001)
    p.add_argument('--sfs', type=float, default=.001)
    p.add_argument('--dt-scale', type=float, default=1.)
    p.add_argument('--dump-every', type=int, default=2000, help='VTP cadence; 0 disables VTP')
    p.add_argument('--snapshot-interval',type=float,default=1.,help='Liquid CSV snapshots for offline angle history; 0 disables')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    if not (math.isfinite(a.end_time) and a.end_time>0 and a.threads>0 and
            math.isfinite(a.sff) and a.sff>0 and math.isfinite(a.sfs) and a.sfs>=0 and
            0<a.dt_scale<=1 and a.dump_every>=0 and
            math.isfinite(a.snapshot_interval) and a.snapshot_interval>=0):
        p.error('Invalid time, threads, strengths, timestep scale or output cadence')
    exe = shutil.which(a.executable)
    if not exe and Path(a.executable).is_file(): exe = str(Path(a.executable).resolve())
    if not exe and not a.dry_run: p.error('Rebuild pf_capillary_validation or supply --executable')
    exe = exe or a.executable
    # Exactly two trajectories, no automated sweeps. Coarse sessile run is a pilot.
    commands = []
    for case, dx in [('droplet', .04), ('sessile', .06)]:
        cmd = [exe, '--case', case, '--dx', str(dx), '--radius', '.6', '--sff', str(a.sff),
               '--sfs', str(a.sfs), '--sound-speed', '10', '--viscosity', '.05',
               '--end-time', str(a.end_time), '--threads', str(a.threads),
               '--dt-scale', str(a.dt_scale), '--snapshot-interval',str(a.snapshot_interval), '--dump-every', str(a.dump_every), '--record-every', '50']
        commands.append((case, cmd))
    if a.dry_run:
        import shlex
        for _, cmd in commands: print(shlex.join(cmd))
        return
    if a.output.exists(): p.error('Choose a new --output directory; existing results are preserved')
    a.output.mkdir(parents=True)
    failures = []
    for case, cmd in commands:
        folder = a.output/case; folder.mkdir()
        (folder/'command.json').write_text(json.dumps(cmd, indent=2)+'\n')
        print(f'Running {case}; log: {folder / "run.log"}', flush=True)
        start = time.monotonic()
        with (folder/'run.log').open('w') as log:
            proc = subprocess.run(cmd, cwd=folder, stdout=log, stderr=subprocess.STDOUT)
        (folder/'execution.json').write_text(json.dumps({'returncode':proc.returncode,
                    'elapsed_seconds':time.monotonic()-start}, indent=2)+'\n')
        if proc.returncode:
            failures.append(case); print(f'{case} failed; partial logs retained', flush=True)
        else:
            print(json.dumps(analyze(folder), indent=2), flush=True)
    # Small diagnostic archive, excluding bulky VTP output. Keep VTP locally for ParaView.
    archive = a.output.resolve().with_suffix('.zip')
    if archive.exists(): p.error('Diagnostic zip already exists; choose another output name')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for case, _ in commands:
            for path in sorted((a.output/case).glob('*')):
                if path.is_file() and path.suffix in ('.csv','.json','.log'):
                    z.write(path, str(Path(a.output.name)/path.relative_to(a.output)))
    print(f'Send back: {archive}', flush=True)
    if failures: raise SystemExit(1)

if __name__ == '__main__':
    main()
