#!/usr/bin/env python3
"""A/B benchmark: SPHinXsys native rebuilds vs cached Verlet candidates.

Runs identical short fixed-sphere bridge initial-value problems, compares the
final particle state and momentum residuals, and archives numerical evidence.
This benchmark does not validate surface tension/contact angle.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import time
import zipfile


def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def run_case(root, exe, mode, args):
    dest = root / mode
    dest.mkdir(parents=True, exist_ok=True)
    cmd = [str(exe), '--case', 'bridge',
           '--dx', str(args.dx), '--radius', str(args.radius),
           '--sphere-radius', str(args.sphere_radius), '--bridge-gap', str(args.bridge_gap),
           '--fixed-spheres', '1', '--sff', str(args.sff), '--sfs', str(args.sfs),
           '--viscosity', str(args.viscosity), '--sound-speed', str(args.sound_speed),
           '--end-time', str(args.end_time), '--steps', str(args.steps),
           '--record-every', '25', '--dump-every', '0', '--threads', str(args.threads),
           '--neighbor-mode', mode, '--neighbor-skin-factor', str(args.skin_factor)]
    start = time.perf_counter()
    result = subprocess.run(cmd, cwd=dest, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - start
    (dest / 'run.log').write_text(result.stdout)
    (dest / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
    if result.returncode:
        raise RuntimeError(f'{mode} failed (exit {result.returncode}); inspect {dest / "run.log"}')
    summary = json.loads((dest / 'pf_run.json').read_text())
    if not summary['completed_time']:
        raise RuntimeError(f'{mode}: run did not reach requested end time')
    history = read_csv(dest / 'pf_history.csv')
    perf = read_csv(dest / 'pf_perf.csv')[-1]
    bridge = read_csv(dest / 'pf_bridge.csv')[-1]
    return {
        'elapsed_seconds': elapsed,
        'solver_seconds': float(perf['wall_seconds']),
        'graph_seconds': float(perf['graph_seconds']),
        'evaluate_seconds': float(perf['evaluate_seconds']),
        'neighbor_builds': int(perf['neighbor_builds']),
        'neighbor_reuses': int(perf['neighbor_reuses']),
        'candidate_pair_evaluations': int(perf['total_pairs']),
        'steps': int(summary['steps']),
        'density_error_final': float(history[-1]['max_density_error']),
        'linear_error_final': float(history[-1]['linear_error']),
        'angular_error_final': float(history[-1]['angular_error']),
        'bridge_force_left_final': float(bridge['left_fx']),
        'bridge_force_right_final': float(bridge['right_fx']),
    }


def compare_state(root):
    a = read_csv(root / 'native' / 'pf_final_liquid.csv')
    b = read_csv(root / 'verlet' / 'pf_final_liquid.csv')
    if len(a) != len(b):
        raise RuntimeError(f'Particle counts mismatch: {len(a)} vs {len(b)}')
    keys = ('x', 'y', 'z', 'rho', 'vx', 'vy', 'vz')
    maxima = {key: 0. for key in keys}
    for ra, rb in zip(a, b):
        for key in keys:
            error = abs(float(ra[key]) - float(rb[key]))
            if not math.isfinite(error):
                raise RuntimeError(f'Nonfinite state comparison for {key}')
            maxima[key] = max(maxima[key], error)
    return {'particles': len(a), 'max_abs_difference': maxima}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('results-pf-neighbor-benchmark'))
    parser.add_argument('--dx', type=float, default=.05)
    parser.add_argument('--radius', type=float, default=.35)
    parser.add_argument('--sphere-radius', type=float, default=.3)
    parser.add_argument('--bridge-gap', type=float, default=.2)
    parser.add_argument('--sff', type=float, default=.001)
    parser.add_argument('--sfs', type=float, default=.002)
    parser.add_argument('--viscosity', type=float, default=.05)
    parser.add_argument('--sound-speed', type=float, default=10.)
    parser.add_argument('--end-time', type=float, default=.15)
    parser.add_argument('--steps', type=int, default=20000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--skin-factor', type=float, default=.25)
    args = parser.parse_args()
    exe = args.executable.resolve()
    if not exe.is_file():
        parser.error(f'Executable not found: {exe}')
    root = args.output.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error(f'Output already exists and is nonempty: {root}; choose another --output')
    root.mkdir(parents=True, exist_ok=True)
    data = {mode: run_case(root, exe, mode, args) for mode in ('native', 'verlet')}
    data['state_difference'] = compare_state(root)
    data['speedup'] = data['native']['solver_seconds'] / data['verlet']['solver_seconds']
    data['wall_speedup'] = data['native']['elapsed_seconds'] / data['verlet']['elapsed_seconds']
    data['matching_steps'] = data['native']['steps'] == data['verlet']['steps']
    maxdiff = data['state_difference']['max_abs_difference']
    data['equivalent_to_tolerance'] = (data['matching_steps'] and
        max(maxdiff.values()) < 1e-7 and
        abs(data['native']['bridge_force_left_final'] - data['verlet']['bridge_force_left_final']) < 1e-7)
    summary = root / 'neighbor_comparison.json'
    summary.write_text(json.dumps(data, indent=2) + '\n')
    archive = root.with_suffix('.zip')
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for item in root.rglob('*'):
            if item.is_file():
                z.write(item, arcname=str(item.relative_to(root.parent)))
    print(json.dumps(data, indent=2))
    print('Results:', archive)
    if not data['equivalent_to_tolerance']:
        raise SystemExit('WARNING: state/force equivalence failed; upload the archive for investigation')


if __name__ == '__main__':
    main()
