#!/usr/bin/env python3
"""Short executable smoke/restart checks. Does not run equilibration or calibration."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


def read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def run(cmd):
    print(' '.join(map(str,cmd)),flush=True)
    subprocess.run(list(map(str,cmd)),check=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',default='build/pf_dem',type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    executable=a.executable.resolve()
    if not executable.is_file():p.error('Executable does not exist')
    if a.output.exists():p.error('Output exists; choose a fresh directory')
    a.output.mkdir(parents=True);out=a.output.resolve()
    runner=[sys.executable,root/'tools/run_pf_dem.py']
    checks=[]
    for name,example in [('two','dem_two_particles.xml'),('levelset','dem_aggregate_levelset.xml')]:
        folder=out/name
        run(runner+[root/'examples'/example,'--executable',executable,'--output',folder,'--end-time','.0002','--steps','20'])
        h=read(folder/'dem_history.csv');b=read(folder/'dem_particles.csv')
        assert json.loads((folder/'dem_run.json').read_text())['completed_time']
        assert max(float(r['linear_error']) for r in h)<1e-11
        assert max(float(r['angular_error']) for r in h)<1e-11
        assert h[0]['liquid_mass']==h[-1]['liquid_mass']
        for file in folder.glob('*.vtp'):ET.parse(file)
        if name=='two':
            dry=next(r for r in b if r['particle']=='1' and r['step']=='0')
            assert all(float(dry['fluid_'+c])==0 for c in ('fx','fy','fz','tx','ty','tz'))
        checks.append(dict(case=name,steps=int(h[-1]['step']),linear_error=float(h[-1]['linear_error']),angular_error=float(h[-1]['angular_error'])))
    # Load a real executable checkpoint and compare its saved endpoint to the
    # resumed start. Exact split-trajectory equivalence is tested in C++.
    xml=out/'restart.xml';xml.write_text('''<pf_dem><fluid dx=".12"/>
    <particle x="0" z="-.45" radius=".24"/><particle x="-.9" radius=".24" vx=".02"/>
    <liquid mode="droplets"><droplet radius=".5"/></liquid>
    <run end_time=".0001" dt_scale=".75" record_every_steps="1" dump_every_steps="0" threads="1"/>
    </pf_dem>''')
    run(runner+[xml,'--executable',executable,'--output',out/'split','--end-time','.0001'])
    run(runner+[xml,'--executable',executable,'--output',out/'resumed','--resume',out/'split/dem_checkpoint.dat','--end-time','.0002'])
    # Exact trajectory restart is covered by the C++ regression; compare the
    # actual executable's loaded start state against its saved endpoint here.
    before=read(out/'split/dem_particles.csv');after=read(out/'resumed/dem_particles.csv')
    for j in ('0','1'):
        x=[r for r in before if r['particle']==j][-1]
        y=next(r for r in after if r['particle']==j)
        for k in ('time','x','y','z','vx','vy','vz','lx','ly','lz','qw','qx','qy','qz'):
            assert x[k]==y[k],(k,x[k],y[k])
    # Restart mismatch must fail before it writes a trajectory.
    xml.write_text(xml.read_text().replace('radius=".24"','radius=".25"',1))
    bad=subprocess.run(runner+[xml,'--executable',executable,'--output',out/'mismatch','--resume',out/'split/dem_checkpoint.dat','--end-time','.0002'],capture_output=True,text=True)
    assert bad.returncode!=0 and 'configuration differs' in bad.stderr
    report=dict(status='short_checks_passed',physics_validated=False,cases=checks,restart_endpoint_exact=True,rejects_changed_geometry=True)
    (out/'checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
