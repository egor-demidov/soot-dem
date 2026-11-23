#!/usr/bin/env python3
"""Minimal actual soot-dem executable startup/restart checks; no equilibration."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import xml.etree.ElementTree as ET


def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',default='build-soot/pf_soot_dem',type=Path)
    p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    exe=a.executable.resolve();out=a.output.resolve();root=Path(__file__).resolve().parents[1]
    if not exe.is_file():p.error('Build pf_soot_dem first')
    if out.exists():p.error('Output exists; choose a fresh directory')
    out.mkdir(parents=True);report=[]
    def run(xml,name,end,resume=None):
        cmd=[sys.executable,str(root/'tools/run_pf_dem.py'),str(xml),'--executable',str(exe),'--output',str(out/name),'--end-time',str(end),'--steps','20']
        if resume:cmd+=['--resume',str(resume)]
        subprocess.run(cmd,check=True)
        h=rows(out/name/'dem_history.csv')
        assert json.loads((out/name/'dem_run.json').read_text())['completed_time']
        assert all(int(r['neck_count'])==2 for r in h)
        assert h[0]['liquid_mass']==h[-1]['liquid_mass']
        assert max(float(r['linear_error']) for r in h)<1e-11
        assert max(float(r['angular_error']) for r in h)<1e-11
        for path in (out/name).glob('necks_*.vtp'):
            piece=ET.parse(path).find('.//Piece');assert piece.get('NumberOfLines')=='2'
        report.append(dict(case=name,steps=int(h[-1]['step']),linear_error=float(h[-1]['linear_error']),angular_error=float(h[-1]['angular_error']),necks=int(h[-1]['neck_count'])))
    for label in ('droplet','levelset'):
        xml=root/f'examples/soot_dem_{label}.xml';run(xml,label,.0002)
        run(xml,label+'_resumed',.0004,out/label/'dem_checkpoint.dat')
        before=rows(out/label/'dem_particles.csv');after=rows(out/(label+'_resumed')/'dem_particles.csv')
        for j in ('0','1','2'):
            x=[r for r in before if r['particle']==j][-1];y=next(r for r in after if r['particle']==j)
            for k in ('time','x','y','z','vx','vy','vz','lx','ly','lz','qw','qx','qy','qz','dem_fx','dem_fy','dem_fz','dem_tx','dem_ty','dem_tz'):
                assert x[k]==y[k],(label,j,k)
    result=dict(status='short_checks_passed',physics_validated=False,cases=report,restart_motion_and_model_loads_exact=True)
    (out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
