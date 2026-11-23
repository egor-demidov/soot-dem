#!/usr/bin/env python3
"""Measure surface tension and apparent contact angle from XML PF-SPH runs."""
import argparse
import json
from pathlib import Path

import numpy as np

from pf_calibration_core import rows, series_statistics


def run_info(folder, case):
    model = json.loads((folder/'resolved_input.json').read_text())['model']
    if model.get('calibration') != case:
        raise ValueError(f'{folder}: expected a {case} calibration run')
    run = json.loads((folder/'dem_run.json').read_text())
    history = rows(folder/'dem_history.csv')
    return model, run, history


def gamma(folders, tail_fraction):
    reports = []
    settings = None
    for folder in folders:
        model, run, history = run_info(folder, 'droplet')
        fluid = model['fluid'][:5]  # sfs is inactive without solid samples.
        if settings is None:
            settings = fluid
        if settings != fluid:
            raise ValueError('Use identical dx, rho0, sound speed, viscosity and sff across radii')
        start, end = float(history[0]['time']), float(history[-1]['time'])
        late = [r for r in history if float(r['time']) >= end-tail_fraction*(end-start)]
        times = [float(r['time']) for r in late]
        if len(times) < 3 or times[0] >= times[-1]:
            raise ValueError(f'{folder}: need at least three late history rows at distinct times')
        min_core = min(int(r['core_particle_count']) for r in late)
        if min_core < 20:
            raise ValueError(f'{folder}: pressure core has {min_core} samples. Increase R/dx')
        p = series_statistics(times, [float(r['virial_pressure']) for r in late])
        radius = series_statistics(times, [float(r['volume_radius']) for r in late])
        reports.append(dict(folder=str(folder), radius=radius, core_pressure=p,
                            core_samples_min=min_core, completed_time=run['completed_time'],
                            tail_max_mach=max(float(r['max_speed']) for r in late)/fluid[2]))
    x = np.array([2/r['radius']['mean'] for r in reports])
    y = np.array([r['core_pressure']['mean'] for r in reports])
    if len(folders) < 3 or len(np.unique(x)) < 3:
        raise ValueError('Provide at least three distinct relaxed droplet radii')
    slope, intercept = np.polyfit(x, y, 1)
    residual = float(np.sqrt(np.mean((y-(slope*x+intercept))**2)))
    failures = []
    for report in reports:
        p = report['core_pressure']
        if (not report['completed_time'] or p['mean']<=0 or report['tail_max_mach']>.02 or
                max(p['drift_over_window'],p['standard_deviation'])>.05*abs(p['mean'])):
            failures.append(f"Pressure or relaxation screen: {report['folder']}")
    if slope<=0 or residual>.1*np.ptp(y) or abs(intercept)>.2*np.mean(abs(y)):
        failures.append('Laplace slope, residual or pressure offset outside screening limits')
    return dict(gamma=float(slope), pressure_offset=float(intercept),
                rms_pressure_residual=residual, screening_failures=failures,
                dx=settings[0], sff=settings[4], droplets=reports)


def angle(folder, frame_count, grid_dx, mesh):
    from pf_interface import analyze_file
    model, run, history = run_info(folder, 'sessile')
    dx, rho0, sound, _, sff, sfs = model['fluid']
    meta = dict(dx=dx, rho0=rho0, case='sessile', wall_z=model['wall']['z'])
    frames = sorted(folder.glob('liquid_[0-9]*.vtp'))[-frame_count:]
    if not frames:
        raise ValueError('No liquid VTP frames. Enable dump_every_steps in the XML')
    by_step = {int(row['step']): row for row in history}
    reports = []
    for frame in frames:
        step = int(frame.stem.split('_')[1])
        row = by_step[step]
        output_mesh = mesh if frame == frames[-1] else None
        try:
            data = analyze_file(frame,meta,grid_dx,output_mesh)
        except ValueError as exc:
            reports.append(dict(step=step,time=float(row['time']),fit_error=str(exc)))
            continue
        reports.append(dict(step=step,time=float(row['time']),mach=float(row['max_speed'])/sound,**data))
    angles = [f['cap_fit']['angle_deg'] for f in reports if f.get('calibration_eligible')]
    angle_range = float(np.ptp(angles)) if len(angles)>1 else None
    passes_screen = (run['completed_time'] and len(angles)>=3 and len(angles)==len(reports)
                     and angle_range<=5 and all(f['mach']<=.02 and
                     f['cap_fit']['angle_cutoff_range_deg']<=8 for f in reports))
    return dict(dx=dx, sff=sff, sfs=sfs, wall_z=meta['wall_z'],
                completed_time=run['completed_time'],
                final_angle_deg=reports[-1].get('cap_fit',{}).get('angle_deg'),
                final_fit_eligible=reports[-1].get('calibration_eligible',False),
                eligible_angle_range_deg=angle_range, passes_stationarity_screen=passes_screen,
                frames=reports)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='case',required=True)
    g = sub.add_parser('gamma')
    g.add_argument('folders',type=Path,nargs='+')
    g.add_argument('--tail-fraction',type=float,default=.25)
    a = sub.add_parser('angle')
    a.add_argument('folder',type=Path)
    a.add_argument('--frames',type=int,default=3)
    a.add_argument('--grid-dx',type=float,default=.25)
    a.add_argument('--mesh',type=Path,help='Write final reconstructed liquid interface as VTP')
    args = parser.parse_args()
    try:
        if args.case=='gamma':
            if not 0<args.tail_fraction<=1 or len(args.folders)<3:
                raise ValueError('Provide at least three runs and a tail fraction in (0, 1]')
            result=gamma(args.folders,args.tail_fraction)
        else:
            if not .05<=args.grid_dx<=1 or args.frames<1:
                raise ValueError('Use at least one frame and grid-dx in [0.05, 1]')
            result=angle(args.folder,args.frames,args.grid_dx,args.mesh)
        print(json.dumps(result,indent=2,allow_nan=False))
    except (KeyError,OSError,ValueError) as exc:
        parser.error(str(exc))


if __name__=='__main__':
    main()
