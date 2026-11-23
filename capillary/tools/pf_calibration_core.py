"""Auditable, conservative postprocessing for PF-SPH calibration runs.

No calibrated value is inferred from a two-step smoke or the paper's s_ff mapping.
Numerical gates are diagnostics, not proofs of continuum accuracy.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import numpy as np


def rows(path):
    with Path(path).open(newline='') as file:
        ans = list(csv.DictReader(file))
    if not ans:
        raise ValueError(f'Empty CSV {path}')
    return ans


def tail(rows_in, fraction=0.2):
    if not 0 < fraction <= 1:
        raise ValueError('Invalid tail fraction')
    final = float(rows_in[-1]['time'])
    if not final > 0:
        raise ValueError('Simulation has no positive elapsed time')
    selected = [r for r in rows_in if float(r['time']) >= final * (1-fraction) - 1e-12]
    if len(selected) < 5:
        raise ValueError('Need >=5 observations in the final time window')
    return selected


def series_statistics(t, y):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(t) < 2 or not np.isfinite(y).all() or not np.isfinite(t).all() or np.ptp(t) <= 0:
        raise ValueError('Insufficient finite and distinct time-series points')
    slope = float(np.polyfit(t-t[0], y, 1)[0])
    return dict(mean=float(y.mean()), standard_deviation=float(y.std(ddof=1)),
                range=float(np.ptp(y)), slope_per_time=slope,
                drift_over_window=abs(slope*float(np.ptp(t))))


def analyze_droplet(folder, min_core=20):
    folder = Path(folder)
    meta = json.loads((folder/'pf_run.json').read_text())
    if meta['case'] != 'droplet':
        raise ValueError('Expected an isolated droplet')
    hist = rows(folder/'pf_history.csv')
    late = tail(hist)
    pressure = [float(r['virial_pressure']) for r in late]
    radii = [float(r['volume_radius']) for r in late]
    core_counts = [int(r.get('core_particle_count', -1)) for r in late]
    if min(core_counts) < min_core:
        raise ValueError(f'Insufficient core particles (minimum {min(core_counts)}; require {min_core}); increase R/dx')
    p_stat = series_statistics([float(r['time']) for r in late], pressure)
    r_stat = series_statistics([float(r['time']) for r in late], radii)
    eos = [float(r['core_eos_pressure']) for r in late]
    eos_stat = series_statistics([float(r['time']) for r in late], eos)
    mach = max(float(r['max_speed']) for r in late)/meta['sound_speed']
    dens = max(float(r['max_density_error']) for r in late)
    relative_drift = p_stat['drift_over_window']/max(abs(p_stat['mean']), 1.e-10)
    relative_std = p_stat['standard_deviation']/max(abs(p_stat['mean']), 1.e-10)
    # A small fitted trend can conceal sustained oscillations. This 5% screen
    # is a conservative numerical gate, not an uncertainty estimate: samples
    # are correlated, and low variation alone cannot establish equilibration.
    checks = dict(completed=bool(meta.get('completed_time')),
                  positive_pressure=math.isfinite(p_stat['mean']) and p_stat['mean'] > 0,
                  pressure_drift=relative_drift <= .05,
                  pressure_fluctuations=relative_std <= .05,
                  mach=mach <= .02, density_error=dens <= .1)
    failures = [name for name, passed in checks.items() if not passed]
    stable = not failures
    return dict(folder=str(folder), dx=meta['dx'], sff=meta['sff'], radius_initial=meta['radius'],
                volume_radius=r_stat['mean'], pressure=p_stat, eos_pressure=eos_stat,
                core_samples_min=min(core_counts), radius_stability=r_stat, tail_mach=mach,
                tail_density_error=dens, pressure_relative_drift=relative_drift,
                pressure_relative_standard_deviation=relative_std,
                stationarity_failures=failures,
                stationary=bool(stable), completed=bool(meta.get('completed_time')),
                note='Virial pressure includes EOS + conservative PF interaction; exterior pressure is zero in this gas-free model.')


def fit_surface_tension(reports):
    if len(reports) < 3:
        raise ValueError('Require >=3 distinct droplet radii to measure pressure versus curvature')
    dx = float(reports[0]['dx'])
    sff = float(reports[0]['sff'])
    if any(not math.isclose(float(a['dx']), dx, abs_tol=1e-12) or
           not math.isclose(float(a['sff']), sff, rel_tol=1e-10) for a in reports):
        raise ValueError('Radii must use identical dx and sff')
    x = 2/np.array([a['volume_radius'] for a in reports], dtype=float)
    y = np.array([a['pressure']['mean'] for a in reports], dtype=float)
    if len(np.unique(np.round(x,12))) != len(x):
        raise ValueError('Duplicate equilibrium droplet radii')
    gamma, offset = np.polyfit(x,y,1)
    fitted = gamma*x+offset
    residual = np.sqrt(np.mean((fitted-y)**2))
    residual_relative = float(residual/max(float(np.ptp(y)),1.e-10))
    offset_relative = float(abs(offset)/max(float(np.mean(abs(y))),1.e-10))
    stable = all(r['stationary'] for r in reports)
    checks = dict(stationary_droplets=stable, positive_gamma=bool(gamma>0),
                  fit_residual=residual_relative<=.1, pressure_offset=offset_relative<=.2)
    failures = [name for name, passed in checks.items() if not passed]
    calibrated = not failures
    return dict(status='provisional_pass' if calibrated else 'not_converged',
                failed_checks=failures,
                dx=dx, sff=sff, gamma=float(gamma), pressure_offset=float(offset),
                rms_fit_residual=float(residual), relative_fit_residual=residual_relative,
                relative_pressure_offset=offset_relative,
                published_model_mapping_comparison=30.96*sff,
                relative_deviation_from_published_mapping=float(gamma/(30.96*sff)-1) if sff else None,
                notes=['Gamma is the fitted slope of total core virial pressure against 2/R.',
                       'The fitted intercept tests for a spurious pressure offset; it is not forced to zero.',
                       'The published 30.96 coefficient is a comparison, not a substituted measurement.',
                       'Passing numerical gates does not establish grid-independent material properties.'],
                droplets=list(reports))


def evaluate_angle_run(folder, grid_dx=.3):
    from pf_interface import analyze_folder
    folder = Path(folder)
    meta = json.loads((folder/'pf_run.json').read_text())
    if meta['case'] != 'sessile':
        raise ValueError('Expected sessile droplet')
    interface = analyze_folder(folder, grid_dx)
    hist = rows(folder/'pf_history.csv')
    late = tail(hist)
    last_fit = interface.get('cap_fit', {})
    history = interface.get('angle_history', {}) or {}
    span = history.get('tail_end',0) - history.get('tail_start',0)
    drift = abs(history.get('tail_angle_slope_deg_per_time',float('inf'))*span)
    valid = (bool(meta.get('completed_time')) and interface.get('calibration_eligible',False) and
             history.get('tail_usable_frames',0)>=3 and
             history.get('tail_angle_range_deg', float('inf')) <=5 and
             drift<=5 and last_fit.get('angle_cutoff_range_deg',float('inf')) <=8 and
             max(float(r['max_density_error']) for r in late)<=.1 and
             max(float(r['max_speed']) for r in late)/meta['sound_speed']<=.02)
    return dict(folder=str(folder), dx=meta['dx'], sff=meta['sff'],sfs=meta['sfs'],
                ratio=meta['sfs']/meta['sff'],angle_deg=last_fit.get('angle_deg'),
                cutoff_range_deg=last_fit.get('angle_cutoff_range_deg'),
                tail_angle_drift_deg=drift,
                tail_angle_range_deg=history.get('tail_angle_range_deg'),
                tail_usable_frames=history.get('tail_usable_frames',0),
                fit_status=last_fit.get('status'),
                volume_error=interface.get('isosurface_volume_relative_error'),
                numerically_stationary=bool(valid),
                note='Angle inferred from volume-constrained, zero-gravity sessile-cap fit. No kinetic or contact-angle hysteresis calibration.')


def select_angle(reports, gamma_report, target_angle, tolerance_deg=5):
    if gamma_report['status']!='provisional_pass':
        raise ValueError('Droplet Laplace-pressure fit must pass before selecting solid interaction')
    valid=[r for r in reports if r['numerically_stationary'] and r['angle_deg'] is not None and
           math.isfinite(r['angle_deg']) and math.isclose(r['dx'], gamma_report['dx'], abs_tol=1e-12) and
           math.isclose(r['sff'], gamma_report['sff'],rel_tol=1e-10)]
    if not valid:
        raise ValueError('No stationary, geometrically resolved sessile-drop fit at matched gamma settings')
    best=min(valid,key=lambda r:abs(r['angle_deg']-target_angle))
    if abs(best['angle_deg']-target_angle)>tolerance_deg:
        raise ValueError(f'Target {target_angle:g} deg not reached: nearest is {best["angle_deg"]:.2f} deg. Refine sfs bracket.')
    return dict(status='provisional_pass',dx=best['dx'],sff=best['sff'],sfs=best['sfs'],
                gamma=gamma_report['gamma'],angle_deg=best['angle_deg'],
                target_angle_deg=target_angle,target_error_deg=best['angle_deg']-target_angle,
                droplet_gamma_report='gamma_report.json',
                note='Provisional calibration at one resolution. Recalibrate gamma and angle at the finer dx, then compare.')


def bridge_summary(folder):
    folder=Path(folder)
    meta=json.loads((folder/'pf_run.json').read_text())
    if meta['case']!='bridge' or not meta.get('fixed_spheres'):
        raise ValueError('Expected fixed-sphere bridge')
    hist=rows(folder/'pf_history.csv')
    bridge=rows(folder/'pf_bridge.csv')
    late=tail(bridge)
    t=[float(r['time']) for r in late]
    fx=series_statistics(t,[float(r['left_fx']) for r in late])
    return dict(dx=meta['dx'],sff=meta['sff'],sfs=meta['sfs'],
                sphere_radius=meta['sphere_radius'],gap=meta['bridge_gap'],liquid_radius=meta['radius'],
                mass=float(hist[-1]['mass']),volume=float(hist[-1]['volume']),
                completed=bool(meta.get('completed_time')),left_force=fx,
                right_force=series_statistics(t,[float(r['right_fx']) for r in late]),
                tail_connected=all(int(r['liquid_connected'])==1 for r in late),
                min_clearance=min(float(r['min_liquid_clearance']) for r in bridge),
                max_density_error=max(float(r['max_density_error']) for r in hist),
                max_mach=max(float(r['max_speed']) for r in hist)/meta['sound_speed'],
                max_linear_residual=max(abs(float(r['linear_error'])) for r in hist),
                max_angular_residual=max(abs(float(r['angular_error'])) for r in hist))


def compare_bridges(coarse,fine,cc,fc,gamma_tol=.1,angle_tol=5.,volume_tol=.04):
    if not(coarse['completed'] and fine['completed'] and coarse['tail_connected'] and fine['tail_connected']):
        raise ValueError('Bridge run incomplete or disconnected')
    if not(coarse['dx']>fine['dx']):
        raise ValueError('Specify coarse then fine resolution')
    for key in ('sphere_radius','gap','liquid_radius'):
        if not math.isclose(coarse[key],fine[key],rel_tol=0,abs_tol=1e-10):
            raise ValueError(f'Mismatched bridge geometry: {key}')
    for r,cal in ((coarse,cc),(fine,fc)):
        if cal['status']!='provisional_pass' or not math.isclose(r['dx'],cal['dx'],abs_tol=1e-12) or not math.isclose(r['sff'],cal['sff'],rel_tol=1e-10) or not math.isclose(r['sfs'],cal['sfs'],rel_tol=1e-10):
            raise ValueError('Bridge run must match its own provisional calibration')
    g_error=abs(cc['gamma']-fc['gamma'])/max(abs(cc['gamma']),abs(fc['gamma']))
    a_error=abs(cc['angle_deg']-fc['angle_deg'])
    v_error=abs(coarse['mass']-fine['mass'])/max(coarse['mass'],fine['mass'])
    if g_error>gamma_tol or a_error>angle_tol or v_error>volume_tol:
        raise ValueError(f'Physical setups are not matched: gamma={g_error:.1%}, angle={a_error:.1f} deg, mass={v_error:.1%}.')
    c=coarse['left_force']['mean'];f=fine['left_force']['mean']
    rel=abs(f-c)/max(abs(c),abs(f),1e-12)
    comparison=dict(status='comparison_only',coarse_force=c,fine_force=f,
                relative_force_difference=rel,
                gamma_relative_difference=g_error,angle_difference_deg=a_error,
                liquid_mass_relative_difference=v_error,
                force_is_converged=bool(rel<=.05),
                note='Two levels alone do not provide an asymptotic convergence order; an independent Young-Laplace reference is needed.')
    angle=.5*(cc['angle_deg']+fc['angle_deg'])
    if 0<angle<180:
        # The reference solves a simple symmetric bridge on the facing hemispheres.
        try:
            from young_laplace import bridge_reference
            ref=bridge_reference(.5*(coarse['volume']+fine['volume']),coarse['gap'],angle,
                                  coarse['sphere_radius'],.5*(cc['gamma']+fc['gamma']))
            comparison['young_laplace_reference']={**ref,
              'relative_fine_force_discrepancy':float((f-ref['attractive_force'])/
                        max(abs(ref['attractive_force']),1e-12))}
        except (ValueError,RuntimeError,FloatingPointError) as exc:
            comparison['young_laplace_reference']={'status':'unresolved','reason':str(exc)}
    else:
        comparison['young_laplace_reference']={'status':'not_supported',
                'reason':'Existing symmetric reference covers contact angles 0<theta<180 degrees'}
    return comparison
