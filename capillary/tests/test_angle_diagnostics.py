import csv
import json
import math
import pathlib
import subprocess
import shutil
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'tools'))
from analyze import exterior_sphere_fit, angle_history, analyze
ROOT=pathlib.Path(__file__).resolve().parents[1]

def cap(theta):
    c=-math.cos(math.radians(theta))
    rows=[]
    for z in np.linspace(max(-1.,-c)+.01,.99,25):
        for phi in np.linspace(0,2*math.pi,32,endpoint=False):
            r=math.sqrt(1-z*z)
            # Include inner layers falsely flagged as surface, as in SPH bands.
            for scale in (1.,.8,.6):
                rows.append(dict(x=scale*r*math.cos(phi),y=scale*r*math.sin(phi),z=c+z,surface=1))
    return rows

class AngleDiagnosticsTests(unittest.TestCase):
    def test_exterior_fit_ignores_interior_surface_band(self):
        for theta in (60,90,120):
            fit=exterior_sphere_fit(cap(theta),.01,'sessile')
            self.assertEqual(fit['status'],'estimated')
            self.assertAlmostEqual(fit['angle_deg'],theta,places=7)
            self.assertLess(fit['fit_rms'],1e-9)
    def test_degenerate_cloud(self):
        rows=[dict(x=i,y=0,z=0) for i in range(20)]
        self.assertEqual(exterior_sphere_fit(rows,.1,'sessile')['status'],'degenerate convex envelope')
    def test_history_tracks_change(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=pathlib.Path(temp)
            with (folder/'liquid_snapshots.csv').open('w') as f:
                f.write('step,time,file\n')
                for i,theta in enumerate((60,65,70,75,80)):
                    name=f'liquid_{i}.csv';f.write(f'{i},{i},{name}\n')
                    with (folder/name).open('w') as g:
                        writer=csv.DictWriter(g,fieldnames=['x','y','z','surface']);writer.writeheader();writer.writerows(cap(theta))
            result=angle_history(folder,dict(case='sessile',dx=.01))
            self.assertAlmostEqual(result['tail_slope_deg_per_time'],5,places=7)
            self.assertFalse(result['equilibrium_established'])
    def test_snapshot_index_is_not_treated_as_particle_snapshot(self):
        fixture=ROOT/'verification/smoke/sessile_dx0.3_theta60_c10_dt1'
        with tempfile.TemporaryDirectory() as temp:
            folder=pathlib.Path(temp)/'case';shutil.copytree(fixture,folder)
            snapshots=sorted(folder.glob('liquid_*.csv'),key=lambda f:int(f.stem.split('_')[-1]))
            last=snapshots[-1]
            meta=json.loads((folder/'run.json').read_text())
            (folder/'liquid_snapshots.csv').write_text(f"step,time,file\n{meta['steps']},{meta['time']},{last.name}\n")
            result=analyze(folder)
            self.assertEqual(result['angle_history']['samples'],1)
            self.assertTrue((folder/'analysis.json').exists())

    def test_focused_runner_has_only_two_cases(self):
        proc=subprocess.run([sys.executable,str(ROOT/'tools/run_contact_angle.py'),'--dry-run'],text=True,capture_output=True,check=True)
        cases=json.loads(proc.stdout)
        self.assertEqual(len(cases),2)
        self.assertEqual([c['name'] for c in cases],['sessile_dx0.06_theta60','sessile_dx0.06_theta120'])
        for case in cases:
            cmd=case['command'];self.assertEqual(cmd[cmd.index('--dump-every')+1],'0')

if __name__=='__main__':unittest.main()
