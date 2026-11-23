#!/usr/bin/env python3
"""Experimental native two-sphere liquid-bridge smoke test and audit.

This verifies native execution and elementary conservation; it does NOT
establish a Young-Laplace solution, a contact angle or a calibrated bridge force.
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


def make_command(a, executable):
    cmd = [str(executable), '--case', 'bridge',
           '--dx', str(a.dx), '--radius', str(a.liquid_radius),
           '--sphere-radius', str(a.sphere_radius), '--sphere-density', str(a.sphere_density),
           '--bridge-gap', str(a.gap), '--fixed-spheres', str(int(a.fixed_spheres)),
           '--contact-stiffness', str(a.contact_stiffness),
           '--sphere-vx', str(a.sphere_vx),
           '--sff', str(a.sff), '--sfs', str(a.sfs),
           '--viscosity', str(a.viscosity), '--sound-speed', str(a.sound_speed),
           '--end-time', str(a.end_time), '--threads', str(a.threads),
           '--record-every', str(a.record_every), '--dump-every', str(a.dump_every),
           '--checkpoint-interval', str(a.end_time), '--steps', str(a.max_steps)]
    if a.resume:
        cmd.extend(['--resume', str(a.resume.resolve())])
    return cmd


def assess(directory, a):
    files = ['pf_history.csv', 'pf_rigid.csv', 'pf_bridge.csv', 'pf_checkpoint.dat',
             'pf_final_liquid.csv', 'pf_run.json']
    for filename in files:
        if not (directory / filename).is_file():
            raise RuntimeError(f'Missing output: {filename}')
    with (directory / 'pf_history.csv').open(newline='') as f:
        history = list(csv.DictReader(f))
    with (directory / 'pf_rigid.csv').open(newline='') as f:
        rigid = list(csv.DictReader(f))
    with (directory / 'pf_bridge.csv').open(newline='') as f:
        bridge = list(csv.DictReader(f))
    summary = json.loads((directory / 'pf_run.json').read_text())
    if len(history) < 2 or len(bridge) < 2 or len(rigid) < 4:
        raise RuntimeError('Too few numerical samples for two-sphere smoke test')
    if not summary['completed_time'] or summary['case'] != 'bridge' or summary['sphere_count'] != 2:
        raise RuntimeError('Native run incomplete or case mismatch')
    final_time = float(bridge[-1]['time'])
    if not math.isclose(final_time, a.end_time, abs_tol=1e-10, rel_tol=0):
        raise RuntimeError(f'Run ended at {final_time}, expected {a.end_time}')
    body_ids = set(row['body_id'] for row in rigid)
    if body_ids != {'0', '1'} or len(rigid) != 2 * len(bridge):
        raise RuntimeError('Did not capture one record for each sphere at every time')
    for row in rigid + bridge + history:
        for name, value in row.items():
            if value is None or not math.isfinite(float(value)) and name not in {'virial_radius', 'virial_pressure'}:
                raise RuntimeError(f'Nonfinite diagnostic {name}={value!r}')
    err_lin = max(abs(float(r['linear_error'])) for r in history)
    err_ang = max(abs(float(r['angular_error'])) for r in history)
    density = max(abs(float(r['max_density_error'])) for r in history)
    clearance = min(float(r['min_liquid_clearance']) for r in bridge)
    neck = min(int(r['neck_particles']) for r in bridge)
    connected = all(int(r['liquid_connected']) == 1 for r in bridge)
    if err_lin > 1e-8 or err_ang > 1e-8:
        raise RuntimeError(f'Conservation diagnostic failed: P={err_lin:g} L={err_ang:g}')
    if clearance < -a.dx * .25:
        raise RuntimeError(f'Liquid penetrated DEM sphere: clearance={clearance:g}')
    if neck == 0 or not connected:
        raise RuntimeError('Bridge is unconnected or has no resolved liquid neck')
    if a.fixed_spheres:
        for row in bridge:
            if not math.isclose(float(row['gap']), a.gap, rel_tol=0, abs_tol=1e-12):
                raise RuntimeError('Fixed-sphere gap changed')
    return {'frames': len(bridge), 'steps': int(history[-1]['step']),
            'linear_error_max': err_lin, 'angular_error_max': err_ang,
            'density_error_max': density, 'clearance_min': clearance,
            'neck_particles_min': neck, 'liquid_connected': connected, 'initial_gap': float(bridge[0]['gap']),
            'final_gap': float(bridge[-1]['gap']),
            'left_force_x_final': float(bridge[-1]['left_fx']),
            'right_force_x_final': float(bridge[-1]['right_fx']),
            'left_direct_pf_x_final': float(bridge[-1]['left_pf_fx']),
            'right_direct_pf_x_final': float(bridge[-1]['right_pf_fx']),
            'physical_validation': 'NOT ESTABLISHED; needs force, angle and resolution benchmarks'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='build/pf_capillary_validation')
    parser.add_argument('--output', type=Path, default=Path('results-pf-bridge-smoke'))
    parser.add_argument('--dx', type=float, default=.1)
    parser.add_argument('--liquid-radius', type=float, default=.35)
    parser.add_argument('--sphere-radius', type=float, default=.3)
    parser.add_argument('--sphere-density', type=float, default=2.)
    parser.add_argument('--gap', type=float, default=.2)
    parser.add_argument('--fixed-spheres', type=int, choices=[0, 1], default=1)
    parser.add_argument('--contact-stiffness', type=float, default=10.)
    parser.add_argument('--sphere-vx', type=float, default=0.)
    parser.add_argument('--sff', type=float, default=.001)
    parser.add_argument('--sfs', type=float, default=.002)
    parser.add_argument('--viscosity', type=float, default=.05)
    parser.add_argument('--sound-speed', type=float, default=10.)
    parser.add_argument('--end-time', type=float, default=.002)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--record-every', type=int, default=1)
    parser.add_argument('--dump-every', type=int, default=0)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    a = parser.parse_args()
    finite = (a.dx, a.liquid_radius, a.sphere_radius, a.sphere_density, a.gap,
              a.contact_stiffness, a.sphere_vx, a.sff, a.sfs, a.viscosity,
              a.sound_speed, a.end_time)
    if (not all(math.isfinite(f) for f in finite) or
        not (a.dx>0 and a.liquid_radius>a.dx and a.sphere_radius>a.dx and
             a.sphere_density>0 and a.gap>0 and a.contact_stiffness>=0 and
             a.sff>=0 and a.sfs>=0 and a.viscosity>=0 and a.sound_speed>0 and
             a.end_time>0 and a.threads>=1 and a.record_every>=1 and
             a.dump_every>=0 and a.max_steps>0)):
        parser.error('Invalid bridge parameters')
    if a.resume and not a.resume.is_file():
        parser.error(f'Checkpoint not found: {a.resume}')
    executable = shutil.which(a.executable)
    if executable is None and Path(a.executable).is_file():
        executable = str(Path(a.executable).resolve())
    if not executable and not a.dry_run:
        parser.error('Build pf_capillary_validation first')
    command = make_command(a, executable or a.executable)
    print(shlex.join(command), flush=True)
    if a.dry_run:
        return
    output = a.output.resolve()
    archive = output.with_suffix('.zip')
    if output.exists() or archive.exists():
        parser.error('Output already exists; use a different --output')
    output.mkdir(parents=True)
    (output / 'command.json').write_text(json.dumps(command, indent=2)+'\n')
    start = time.monotonic()
    with (output / 'run.log').open('w') as log:
        result = subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT)
    report = {'returncode': result.returncode, 'wall_seconds': time.monotonic()-start}
    if result.returncode == 0:
        try:
            report.update(assess(output, a))
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            report['error'] = str(exc)
    (output / 'check.json').write_text(json.dumps(report, indent=2)+'\n')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in sorted(output.iterdir()):
            if file.is_file() and file.suffix in {'.csv', '.json', '.log', '.dat'}:
                zipf.write(file, str(Path(output.name) / file.name))
    print(f'Archive: {archive}; check: {report}')
    if result.returncode != 0 or 'error' in report:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
