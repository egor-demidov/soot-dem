"""XML/aggregate preparation checks, without running fluid trajectories."""
from pathlib import Path
import sys
import tempfile
import unittest
import subprocess
import json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_pf_dem import parse_xml, prepared_text
ROOT=Path(__file__).resolve().parents[1]


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.tmp_path=Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_examples_and_initial_dry_particle(self):
        c=parse_xml(ROOT/'examples/dem_two_particles.xml');m=c['model']
        assert len(m['bodies'])==2
        a=m['bodies'][1];d=m['drops'][0]
        distance=sum((a['center'][k]-d['center'][k])**2 for k in range(3))**.5
        assert distance-a['radius']-d['radius']>4*m['fluid'][0]
        assert a['velocity'][0]>0
        assert prepared_text(c).startswith('PFDEM_INPUT 1\n'+c['signature'])
        level=parse_xml(ROOT/'examples/dem_aggregate_levelset.xml')
        assert len(level['model']['bodies'])==3 and level['model']['rm']==.25


    def test_vtk_aggregate_positions_and_radii(self):
        tmp_path=self.tmp_path
        (tmp_path/'aggregate.vtk').write_text('''# vtk DataFile Version 3.0
    aggregate
    ASCII
    DATASET POLYDATA
    POINTS 2 double
    -1 0 0 1 0 0
    POINT_DATA 2
    SCALARS radius double 1
    LOOKUP_TABLE default
    .2 .3
    VECTORS velocity double
    9 9 9 8 8 8
    ''')
        xml=tmp_path/'run.xml';xml.write_text('''<pf_dem><fluid dx=".1"/>
        <particles file="aggregate.vtk" density="2" length_scale="2"/>
        <liquid mode="droplets"><droplet radius="3"/></liquid></pf_dem>''')
        c=parse_xml(xml);a,b=c['model']['bodies']
        assert a['center']==[-2,0,0] and b['radius']==.6
        assert a['velocity']==b['velocity']==[0,0,0]


    def test_particle_csv_velocities_ignored_by_default(self):
        tmp_path=self.tmp_path
        (tmp_path/'aggregate.csv').write_text('x,y,z,radius,vx,vy,vz\n0,0,0,.2,5,6,7\n')
        path=tmp_path/'run.xml';text='''<pf_dem><fluid dx=".05"/>
        <particles file="aggregate.csv"/><liquid mode="droplets"><droplet radius=".5"/></liquid>
        <run end_time="1"/></pf_dem>''';path.write_text(text)
        first=parse_xml(path)
        assert first['model']['bodies'][0]['velocity']==[0,0,0]
        path.write_text(text.replace('end_time="1"','end_time="3"'))
        assert parse_xml(path)['signature']==first['signature']
        path.write_text(text.replace('file="aggregate.csv"','file="aggregate.csv" ignore_velocities="false"'))
        second=parse_xml(path)
        assert second['model']['bodies'][0]['velocity']==[5,6,7]
        assert second['signature']!=first['signature']


    def test_rejects_bad_geometry_and_typos(self):
        tmp_path=self.tmp_path
        p=tmp_path/'bad.xml'
        p.write_text('<pf_dem><fluid dx="nan"/></pf_dem>')
        with self.assertRaisesRegex(ValueError,'dx'):parse_xml(p)
        p.write_text('<pf_dem><fluid dx=".1" typo="1"/></pf_dem>')
        with self.assertRaisesRegex(ValueError,'Unknown'):parse_xml(p)
        p.write_text('<pf_dem><fluid/><particle radius=".2"/><liquid mode="levelset" threshold="2"/></pf_dem>')
        with self.assertRaisesRegex(ValueError,'threshold'):parse_xml(p)

    def test_relative_executable_survives_output_directory_change(self):
        executable=self.tmp_path/'build/pf_dem'
        executable.parent.mkdir()
        executable.write_text('#!/bin/sh\nexit 0\n')
        executable.chmod(0o755)
        command=[sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(ROOT/'examples/dem_two_particles.xml'),
                 '--executable','build/pf_dem','--output','results']
        run=subprocess.run(command,cwd=self.tmp_path,capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr)
        saved=json.loads((self.tmp_path/'results/command.json').read_text())
        self.assertEqual(saved[0],str(executable.resolve()))

if __name__=='__main__':unittest.main()
