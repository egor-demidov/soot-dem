#!/usr/bin/env python3
"""Recover an existing dx=.04 sessile wetting CSV, run a SHORT verification,
then preserve a full-state checkpoint. No validation/equilibrium claim.
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


def read_snapshot_index(path):
    index = path.parent / 'pf_snapshots.csv'
    if not index.is_file():
        raise ValueError(f'Missing {index}')
    with index.open(newline='') as f:
        rows = list(csv.DictReader(f))
    match = [r for r in rows if r['file'] == path.name]
    if len(match) != 1:
        raise ValueError('Snapshot must occur exactly once in pf_snapshots.csv')
    return int(match[0]['step']), float(match[0]['time'])


def make_command(a, executable, snapshot):
    return [executable, '--case', 'sessile', '--dx', str(a.dx), '--radius', str(a.radius),
            '--sff', str(a.sff), '--sfs', str(a.sfs), '--viscosity', str(a.viscosity),
            '--sound-speed', str(a.sound_speed), '--dt-scale', str(a.dt_scale),
            '--end-time', str(a.end_time), '--threads', str(a.threads),
            '--record-every', '50', '--snapshot-interval', '0',
            '--dump-every', '0', '--checkpoint-interval', '.02',
            '--resume-snapshot', str(snapshot)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot', type=Path, required=True, help='pf_snapshot_23377.csv beside pf_snapshots.csv')
    p.add_argument('--executable', default='build/pf_capillary_validation')
    p.add_argument('--output', type=Path, default=Path('results-pf-dx004-restart-smoke'))
    p.add_argument('--end-time', type=float, default=7.05,
                   help='Absolute target simulation time, NOT time after recovery')
    p.add_argument('--threads', type=int, default=6)
    p.add_argument('--dx', type=float, default=.04)
    p.add_argument('--radius', type=float, default=.6)
    p.add_argument('--sff', type=float, default=.001)
    p.add_argument('--sfs', type=float, default=.0035)
    p.add_argument('--viscosity', type=float, default=.05)
    p.add_argument('--sound-speed', type=float, default=20.,
                   help='Inferred from saved dt/history; verify if the original command differed')
    p.add_argument('--dt-scale', type=float, default=1.,
                   help='Inferred from saved dt/history; verify if the original command differed')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    snapshot = a.snapshot.resolve()
    if not snapshot.is_file():
        p.error('Snapshot CSV does not exist')
    try:
        step, starting_time = read_snapshot_index(snapshot)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        p.error(str(exc))
    if not (math.isfinite(a.end_time) and a.end_time > starting_time and
            a.threads > 0 and a.dx > 0 and a.radius > a.dx and
            all(math.isfinite(v) and v >= 0 for v in
                (a.sff,a.sfs,a.viscosity)) and a.sound_speed > 0 and 0 < a.dt_scale <= 1):
        p.error('Invalid simulation options or end-time before snapshot')
    executable = shutil.which(a.executable)
    if not executable and Path(a.executable).is_file():
        executable = str(Path(a.executable).resolve())
    if not executable and not a.dry_run:
        p.error('Executable not found: build pf_capillary_validation')
    cmd = make_command(a, executable or a.executable, snapshot)
    print(f'Recovering snapshot at t={starting_time:.12f}, step {step}')
    print(shlex.join(cmd), flush=True)
    if a.dry_run:
        return
    target = a.output.resolve()
    archive = target.with_suffix('.zip')
    if target.exists() or archive.exists():
        p.error('Output folder/archive already exists; choose a new --output')
    target.mkdir(parents=True)
    (target / 'command.json').write_text(json.dumps(cmd, indent=2)+'\n')
    meta = dict(snapshot=str(snapshot), source_step=step, source_time=starting_time,
                recovery_ledger='Rebased at snapshot; historical wall impulse unavailable',
                original_sound_speed_and_dt_scale='Inferred from observed acoustic step; not file-backed',
                physics_status='NOT_VALIDATED')
    (target / 'recovery.json').write_text(json.dumps(meta, indent=2)+'\n')
    start = time.monotonic()
    with (target / 'run.log').open('w') as log:
        proc = subprocess.run(cmd, cwd=target, stdout=log, stderr=subprocess.STDOUT)
    (target / 'execution.json').write_text(json.dumps(dict(returncode=proc.returncode,
        wall_seconds=time.monotonic()-start),indent=2)+'\n')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for file in sorted(target.rglob('*')):
            if file.is_file() and file.suffix in ('.csv','.json','.log','.dat'):
                z.write(file, str(Path(target.name) / file.relative_to(target)))
    print(f'Exit code={proc.returncode}. Results (including checkpoint): {archive}')
    if proc.returncode:
        raise SystemExit(proc.returncode)


if __name__=='__main__':
    main()
