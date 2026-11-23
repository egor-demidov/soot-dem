"""VTP reader/writer and series metadata checks; no fluid simulations."""
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from export_dem_surface import read_liquid,write_surface,export
from pf_interface import density_grid

class SurfaceIOTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def points(self,step,xyz):
        file=self.root/f'liquid_{step:010d}.vtp'
        r=ET.Element('VTKFile',type='PolyData');poly=ET.SubElement(r,'PolyData')
        piece=ET.SubElement(poly,'Piece',NumberOfPoints=str(len(xyz)))
        a=ET.SubElement(ET.SubElement(piece,'Points'),'DataArray',type='Float64',format='ascii',NumberOfComponents='3')
        a.text=' '.join(str(v) for row in xyz for v in row)
        ET.ElementTree(r).write(file);return file
    def test_point_cloud_roundtrip_and_compact_support(self):
        points=np.array([[0,0,0],[.1,0,0],[0,.1,0]])
        p,time=read_liquid(self.points(0,points))
        np.testing.assert_array_equal(p,points);self.assertIsNone(time)
        grid,origin=density_grid(p,np.full(3,.001),.4,.05)
        self.assertTrue(np.isfinite(grid).all());self.assertGreater(grid.max(),0)
        # Integrated compact kernel mass tends to sample mass with fine grid.
        self.assertAlmostEqual(float(grid.sum())*.05**3,.003,delta=2e-6)
    def test_surface_normals_and_connectivity(self):
        xyz=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        f=self.root/'surface.vtp';write_surface(f,xyz,np.array([[0,1,2]]),np.tile([0.,0,1],(3,1)))
        p=ET.parse(f).find('.//Piece');self.assertEqual(p.get('NumberOfPolys'),'1')
        self.assertEqual(p.find('PointData').get('Normals'),'Normals')
        self.assertEqual(p.find("Polys/DataArray[@Name='connectivity']").text.split(),['0','1','2'])
    def test_series_times_stride_and_output_guard(self):
        (self.root/'resolved_input.json').write_text(json.dumps({'model':{'fluid':[.1,1]}}))
        (self.root/'dem_history.csv').write_text('step,time\n0,0\n2,.001\n4,.003\n')
        xyz=[[x,y,z] for x in (-.1,0,.1) for y in (-.1,0,.1) for z in (-.1,0,.1)]
        for step in (0,2,4):self.points(step,xyz)
        # Isolate series orchestration from the established third-party mesher.
        # This environment lacks scikit-image; this is not a mesher validation.
        mod=types.ModuleType('skimage.measure')
        def mesh(*args,**kwargs):
            return np.array([[0.,0,0],[1,0,0],[0,1,0]]),np.array([[0,1,2]]),np.tile([0.,0,1],(3,1)),None
        mod.marching_cubes=mesh
        with patch.dict(sys.modules,{'skimage':types.ModuleType('skimage'),'skimage.measure':mod}):
            result=export(self.root,self.root/'surfaces',every=2)
            self.assertEqual([r['step'] for r in result['frames']],[0,4])
            self.assertEqual(result['animation_time'],'simulation_time')
            self.assertTrue(all(r['resolved'] and r['vertices']==3 for r in result['frames']))
            datasets=ET.parse(self.root/'surfaces/liquid_surface.pvd').findall('.//DataSet')
            self.assertEqual([float(r.get('timestep')) for r in datasets],[0,.003])
            with self.assertRaisesRegex(ValueError,'exists'):export(self.root,self.root/'surfaces')

if __name__=='__main__':unittest.main()
