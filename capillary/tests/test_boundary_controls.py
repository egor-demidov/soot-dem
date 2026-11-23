from pathlib import Path
import sys,tempfile,unittest,subprocess,json,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_boundary_controls import validate_base,make_control
from run_pf_dem import parse_xml
ROOT=Path(__file__).resolve().parents[1]
class BoundaryControlsTests(unittest.TestCase):
 def test_control_changes_only_solid_strength_and_checkpoint_cadence(self):
  root,base=validate_base(ROOT/'examples/boundary_fixed_sphere.xml')
  for strength in (0.,.006):
   with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'case.xml';ET.ElementTree(make_control(root,strength)).write(p)
    c=parse_xml(p);expected=dict(base['model']);expected['fluid']=base['model']['fluid'][:];expected['fluid'][5]=strength
    self.assertEqual(c['model'],expected)
    self.assertEqual(c['run']['end_time'],4e-10)
    self.assertTrue(c['model']['bodies'][0]['fixed'])
 def test_free_sphere_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'case.xml';p.write_text((ROOT/'examples/boundary_fixed_sphere.xml').read_text().replace('fixed="true"','fixed="false"'))
   with self.assertRaisesRegex(ValueError,'fixed'):validate_base(p)
 def test_checkpoint_audit_simple_reference(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'cp.dat';out=Path(td)/'audit.json'
   # One sphere, two liquid samples, one rigid interior sample. No mechanics.
   rows=['PFDEM 2','"'+'a'*64+'"','.1 1000 10 .001 .002 .012','0 0 0 .01','2 3 1 1 0 0']+['0 0 0']*6
   rows+=['0 0 0 0 0 0 0 0 0 1 0 0 0 .2 1 1',
          '.25 0 0 0 0 0 800 .001 0','.3 0 0 0 0 0 1100 .001 0',
          '.1 0 0 0 0 0 1000 .001 1','2 0 .1 0 0']
   p.write_text('\n'.join(rows)+'\n')
   result=subprocess.run([sys.executable,str(ROOT/'tools/audit_boundary_checkpoint.py'),str(p),'--output',str(out)],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   data=json.loads(out.read_text());self.assertEqual(data['inside_fluid_count'],0)
   self.assertAlmostEqual(data['probes'][0]['rho_ratio'],.8)
   self.assertEqual(data['probes'][0]['components']['solid']['count'],1)
   self.assertGreater(data['probes'][0]['components']['solid']['pair_force']['outward_force_N'],0) # q=.375: repulsive fs
   # Both extrema can legitimately lie outside every boundary kernel support.
   out.unlink();p.write_text(('\n'.join(rows)+'\n').replace('.25 0 0 0 0 0 800','.8 0 0 0 0 0 800').replace('.3 0 0 0 0 0 1100','.9 0 0 0 0 0 1100'))
   result=subprocess.run([sys.executable,str(ROOT/'tools/audit_boundary_checkpoint.py'),str(p),'--output',str(out)],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   data=json.loads(out.read_text())
   for probe in data['probes']:
    self.assertIsNone(probe['nearest_boundary_sample_distance_dx'])
    self.assertEqual(probe['components']['solid']['count'],0)
    self.assertEqual(probe['components']['solid']['pressure']['force_N'],[0.,0.,0.])
   # A sample at a sphere center has no defined surface normal: JSON null.
   out.unlink();p.write_text(('\n'.join(rows)+'\n').replace('.25 0 0 0 0 0 800','0 0 0 0 0 0 800'))
   result=subprocess.run([sys.executable,str(ROOT/'tools/audit_boundary_checkpoint.py'),str(p),'--output',str(out)],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   self.assertIsNone(json.loads(out.read_text())['probes'][0]['relative_outward_velocity_m_s'])

if __name__=='__main__':unittest.main()
