"""Physical-parameter matching and numerical stationarity regression checks."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from pf_calibration_core import (series_statistics,analyze_droplet,fit_surface_tension,
    select_angle,compare_bridges,bridge_summary)


def droplet(radius,dx=.05,sff=.001,gamma=.03,offset=.001,stationary=True):
    return dict(dx=dx,sff=sff,volume_radius=radius,
                pressure=dict(mean=2*gamma/radius+offset),stationary=stationary)


def test_laplace_fit_unbiased_three_radii():
    data=[droplet(r,gamma=.0375,offset=.0008) for r in [.65,.8,.95]]
    fit=fit_surface_tension(data)
    assert fit['status']=='provisional_pass'
    assert abs(fit['gamma']-.0375)<1e-12
    assert abs(fit['pressure_offset']-.0008)<1e-12


def test_laplace_fit_rejects_short_unstable_or_mismatched():
    with pytest.raises(ValueError,match='>=3'):
        fit_surface_tension([droplet(.65),droplet(.8)])
    with pytest.raises(ValueError,match='identical'):
        fit_surface_tension([droplet(.65),droplet(.8,dx=.025),droplet(.95)])
    bad=[droplet(.65),droplet(.8),droplet(.95,stationary=False)]
    assert fit_surface_tension(bad)['status']=='not_converged'


def test_droplet_core_requires_enough_particles(tmp_path):
    meta={'case':'droplet','dx':.05,'sff':.001,'radius':.8,'time':10,
          'sound_speed':10,'completed_time':True}
    (tmp_path/'pf_run.json').write_text(json.dumps(meta))
    with (tmp_path/'pf_history.csv').open('w',newline='') as f:
        fields=['time','virial_pressure','volume_radius','core_particle_count',
                'core_eos_pressure','max_speed','max_density_error']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for t in np.linspace(0,10,21):
            w.writerow(dict(time=t,virial_pressure=.075,volume_radius=.8,
                           core_particle_count=15,core_eos_pressure=.02,max_speed=.01,max_density_error=.03))
    with pytest.raises(ValueError,match='Insufficient core'):
        analyze_droplet(tmp_path)
    text=(tmp_path/'pf_history.csv').read_text().replace(',15,',',30,')
    (tmp_path/'pf_history.csv').write_text(text)
    got=analyze_droplet(tmp_path)
    assert got['stationary']
    assert abs(got['pressure']['mean']-.075)<1e-12


def test_stationarity_rejects_pressure_drift(tmp_path):
    meta={'case':'droplet','dx':.05,'sff':.001,'radius':.8,'time':10,
          'sound_speed':10,'completed_time':True}
    (tmp_path/'pf_run.json').write_text(json.dumps(meta))
    with (tmp_path/'pf_history.csv').open('w',newline='') as f:
        fields=['time','virial_pressure','volume_radius','core_particle_count',
                'core_eos_pressure','max_speed','max_density_error']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for t in np.linspace(0,10,41):
            w.writerow(dict(time=t,virial_pressure=.05+0.03*t,volume_radius=.8,
                           core_particle_count=100,core_eos_pressure=.01,max_speed=.01,max_density_error=.03))
    assert not analyze_droplet(tmp_path)['stationary']


def test_angle_selection_doesnt_silently_pick_unmatched():
    gamma=dict(status='provisional_pass',dx=.05,sff=.001,gamma=.031)
    angles=[dict(dx=.05,sff=.001,sfs=s,angle_deg=t,numerically_stationary=True) for s,t in
            [(.002,110),(.003,87),(.004,63)]]
    result=select_angle(angles,gamma,60)
    assert math.isclose(result['sfs'],.004)
    assert result['status']=='provisional_pass'
    with pytest.raises(ValueError,match='Target'):
        select_angle(angles,gamma,35)
    angles[2]['numerically_stationary']=False
    with pytest.raises(ValueError,match='Target'):
        select_angle(angles,gamma,60)


def bridge(dx, mass=.17, force=-.05):
    return dict(completed=True,tail_connected=True,dx=dx,sphere_radius=.3,gap=.2,
                liquid_radius=.35,mass=mass,volume=mass/1.03,left_force=dict(mean=force),sff=.001,sfs=.002)


def cal(dx, gamma=.03,angle=60):
    return dict(status='provisional_pass',dx=dx,sff=.001,sfs=.002,gamma=gamma,angle_deg=angle)


def test_comparison_checks_physical_matching():
    coarse=bridge(.05);fine=bridge(.025,force=-.048)
    a=compare_bridges(coarse,fine,cal(.05),cal(.025))
    assert a['force_is_converged']
    assert a['relative_force_difference']<.05
    with pytest.raises(ValueError,match='Physical setups'):
        compare_bridges(coarse,fine,cal(.05,angle=60),cal(.025,angle=80))
    with pytest.raises(ValueError,match='Physical setups'):
        compare_bridges(coarse,bridge(.025,mass=.13),cal(.05),cal(.025))
    with pytest.raises(ValueError,match='match its own'):
        compare_bridges(coarse,bridge(.025),cal(.05),cal(.05))


def test_calibration_dry_run_commands():
    cmd=[sys.executable,str(ROOT/'tools/pf_calibrate.py'),'droplets', '--dx','.05',
         '--radii','.65','.8','.95','--output','unused','--dry-run','--executable','nonexisting']
    p=subprocess.run(cmd,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    assert p.stdout.count('--case droplet')==3
    assert '--snapshot-interval 1.0' in p.stdout
    bad=subprocess.run(cmd[:-2]+['--radii','.4'],capture_output=True,text=True)
    assert bad.returncode!=0


def test_stationarity_rejects_zero_trend_oscillations(tmp_path):
    meta={'case':'droplet','dx':.05,'sff':.001,'radius':.8,'time':10,
          'sound_speed':10,'completed_time':True}
    (tmp_path/'pf_run.json').write_text(json.dumps(meta))
    with (tmp_path/'pf_history.csv').open('w',newline='') as f:
        fields=['time','virial_pressure','volume_radius','core_particle_count',
                'core_eos_pressure','max_speed','max_density_error']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        # Final window has symmetric oscillations, with zero fitted trend.
        for k,t in enumerate(np.linspace(0,10,41)):
            w.writerow(dict(time=t,virial_pressure=.1+.03*(-1)**k,volume_radius=.8,
                           core_particle_count=100,core_eos_pressure=4,max_speed=.001,max_density_error=.04))
    report=analyze_droplet(tmp_path)
    assert report['pressure_relative_drift']<1e-12
    assert report['pressure_relative_standard_deviation']>.25
    assert report['stationarity_failures']==['pressure_fluctuations']
    assert not report['stationary']


def test_fit_explains_unsteady_cases_and_nonzero_offset():
    data=[droplet(r,offset=-.02) for r in [.65,.8,.95]]
    data[-1]['stationary']=False
    result=fit_surface_tension(data)
    assert result['failed_checks']==['stationary_droplets','pressure_offset']
    assert result['status']=='not_converged'
    assert abs(result['gamma']-.03)<1e-12
