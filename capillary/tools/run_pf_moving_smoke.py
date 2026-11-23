#!/usr/bin/env python3
"""Short *experimental* moving DEM sphere + PF liquid native SPHinXsys check.

Does not validate wettability/contact angle or a formed capillary bridge.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import zipfile


def command(a, exe):
    return [exe, '--case', 'immersed', '--dx', str(a.dx),
            '--radius', str(a.radius), '--sphere-radius', str(a.sphere_radius),
            '--sphere-density', str(a.sphere_density), '--sphere-vx', str(a.sphere_vx),
            '--sphere-omega-z', str(a.sphere_omega_z),
            '--sff', str(a.sff), '--sfs', str(a.sfs),
            '--viscosity', str(a.viscosity), '--sound-speed', str(a.sound_speed),
            '--end-time', str(a.end_time), '--threads', str(a.threads),
            '--record-every', '1', '--dump-every', '0',
            '--checkpoint-interval', str(a.end_time)]


def check_outputs(folder, end_time, sphere_vx):
    rigid = folder / 'pf_rigid.csv'
    history = folder / 'pf_history.csv'
    checkpoint = folder / 'pf_checkpoint.dat'
    for f in (rigid, history, checkpoint):
        if not f.is_file():
            raise RuntimeError(f'Missing output: {f}')
    with rigid.open(newline='') as f:
        bodies = list(csv.DictReader(f))
    with history.open(newline='') as f:
        hist = list(csv.DictReader(f))
    if len(bodies) < 2 or len(hist) < 2:
        raise RuntimeError('Too few recorded states for moving DEM validation')
    latest = hist[-1]
    if not math.isclose(float(latest['time']), end_time, rel_tol=0, abs_tol=1e-10):
        raise RuntimeError(f'Native simulation did not reach end time: {latest["time"]}')
    fields = ('linear_error', 'angular_error', 'max_density_error')
    residuals = {key: float(latest[key]) for key in fields}
    if not all(math.isfinite(v) for v in residuals.values()):
        raise RuntimeError('Nonfinite conserved-momentum or density diagnostics')
    for row in bodies:
        if not all(math.isfinite(float(value)) for key, value in row.items() if key != 'step'):
            raise RuntimeError('Nonfinite rigid-body diagnostic')
    initial_x = float(bodies[0]['cx'])
    final_x = float(bodies[-1]['cx'])
    if sphere_vx != 0 and abs(final_x - initial_x) < 1e-12:
        raise RuntimeError('Rigid sphere never translated; check force integration')
    return dict(residuals, initial_x=initial_x, final_x=final_x, frames=len(bodies))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='build/pf_capillary_validation')
    parser.add_argument('--output', type=Path, default=Path('results-pf-moving-smoke'))
    parser.add_argument('--dx', type=float, default=.12)
    parser.add_argument('--radius', type=float, default=.8)
    parser.add_argument('--sphere-radius', type=float, default=.24)
    parser.add_argument('--sphere-density', type=float, default=2.)
    parser.add_argument('--sphere-vx', type=float, default=.05)
    parser.add_argument('--sphere-omega-z', type=float, default=.5)
    parser.add_argument('--sff', type=float, default=.001)
    parser.add_argument('--sfs', type=float, default=.002)
    parser.add_argument('--viscosity', type=float, default=.05)
    parser.add_argument('--sound-speed', type=float, default=10.)
    parser.add_argument('--end-time', type=float, default=.002)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--dry-run', action='store_true')
    a = parser.parse_args()
    if not (a.dx > 0 and a.radius > a.sphere_radius + 2*a.dx and
            a.sphere_density > 0 and a.threads >= 1 and a.end_time > 0 and
            a.sound_speed > 0 and all(math.isfinite(x) for x in (
            a.dx, a.radius, a.sphere_radius, a.sphere_density, a.sphere_vx,
            a.sphere_omega_z, a.sff, a.sfs, a.viscosity, a.sound_speed, a.end_time)) and
            a.sff >= 0 and a.sfs >= 0 and a.viscosity >= 0):
        parser.error('Invalid parameters')
    executable = shutil.which(a.executable)
    if executable is None and Path(a.executable).is_file():
        executable = str(Path(a.executable).resolve())
    if not executable and not a.dry_run:
        parser.error('Native executable missing: build pf_capillary_validation')
    cmd = command(a, executable or a.executable)
    print(shlex.join(cmd), flush=True)
    if a.dry_run:
        return
    folder = a.output.resolve()
    archive = folder.with_suffix('.zip')
    if folder.exists() or archive.exists():
        parser.error('Output already exists; use a new --output')
    folder.mkdir(parents=True)
    (folder / 'command.json').write_text(json.dumps(cmd, indent=2)+'\n')
    begin = time.monotonic()
    with (folder / 'run.log').open('w') as log:
        result = subprocess.run(cmd, cwd=folder, stdout=log, stderr=subprocess.STDOUT)
    report = dict(returncode=result.returncode, wall_seconds=time.monotonic()-begin,
                  validation_status='NOT validated physically; short numerical smoke only')
    if result.returncode == 0:
        try:
            report.update(check_outputs(folder,a.end_time,a.sphere_vx))
        except (ValueError, KeyError, RuntimeError) as exc:
            report['error'] = str(exc)
    (folder / 'check.json').write_text(json.dumps(report,indent=2)+'\n')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in sorted(folder.iterdir()):
            if f.is_file() and f.suffix in {'.csv','.json','.log','.dat'}:
                z.write(f, str(Path(folder.name) / f.name))
    print(f'Results: {archive}. Diagnostics: {report}')
    if result.returncode:
        raise SystemExit(result.returncode)
    if 'error' in report:
        raise SystemExit(report['error'])


if __name__ == '__main__':
    main()
