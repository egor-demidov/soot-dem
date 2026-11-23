from pathlib import Path
import unittest,sys,tempfile,subprocess,os
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_pf_dem import parse_xml,prepared_text
ROOT=Path(__file__).resolve().parents[1]
class DensityCheckOptionTests(unittest.TestCase):
 def test_xml_flag_is_runtime_only(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'input.xml';base=(ROOT/'examples/boundary_fixed_sphere.xml').read_text()
   p.write_text(base);a=parse_xml(p)
   p.write_text(base.replace('<run ','<run density_check="false" '));b=parse_xml(p)
   self.assertTrue(a['run']['density_check']);self.assertFalse(b['run']['density_check'])
   self.assertEqual(a['signature'],b['signature']);self.assertEqual(prepared_text(a),prepared_text(b))
 def test_cli_flag_forwarded_without_new_input_version(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td);out=td/'out'
   r=subprocess.run([sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(ROOT/'examples/boundary_fixed_sphere.xml'),'--output',str(out),'--prepare-only','--disable-density-check'],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr)
   import json
   self.assertIn('--disable-density-check',json.loads((out/'command.json').read_text()))
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for tiny native stop/continue checks')
 def test_native_density_disabled_but_mach_retained(self):
  exe=os.environ['PF_DEM_TEST_EXE']
  with tempfile.TemporaryDirectory() as td:
   td=Path(td);xml=td/'in.xml'
   xml.write_text('<pf_dem><fluid dx=".1" rho0="1000" sound_speed="10"/><particle radius=".2" density="2" fixed="true"/><liquid mode="droplets"><droplet radius=".5"/></liquid></pf_dem>')
   config=parse_xml(xml);(td/'prepared.dat').write_text(prepared_text(config))
   rows=['PFDEM 2','"'+config['signature']+'"','.1 1000 10 .05 .001 .003','100 0 0 .01','2 3 1 1 0 0']+['0 0 0']*6
   rows+=['0 0 0 0 0 0 0 0 0 1 0 0 0 .2 1 1','.25 0 0 0 0 0 690 .001 0','.3 0 0 0 0 0 1000 .001 0','.1 0 0 0 0 0 1000 .001 1','2 0 .1 0 0','"none"']
   cp=td/'cp.dat';cp.write_text('\n'.join(rows)+'\n')
   def run(name,flags):
    cwd=td/name;cwd.mkdir()
    return subprocess.run([exe,'--input',str(td/'prepared.dat'),'--resume',str(cp),'--end-time','1e-12','--steps','1']+flags,cwd=cwd,capture_output=True,text=True)
   self.assertEqual(run('enabled',[]).returncode,1)
   a=run('disabled',['--disable-density-check']);self.assertEqual(a.returncode,0,a.stderr)
   # Startup speed exceeds Mach 0.2 even with density check disabled.
   cp.write_text(('\n'.join(rows)+'\n').replace('.25 0 0 0 0 0 690','.25 0 0 3 0 0 690'))
   a=run('mach',['--disable-density-check']);self.assertEqual(a.returncode,1);self.assertIn('Mach=',a.stderr)
if __name__=='__main__':unittest.main()
