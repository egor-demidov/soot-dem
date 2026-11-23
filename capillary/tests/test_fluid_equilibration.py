"""XML, event timing and checkpoint continuation for the initial fluid hold."""
from pathlib import Path
import csv,json,math,os,subprocess,sys,tempfile,unittest
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_pf_dem import parse_xml,prepared_text
ROOT=Path(__file__).resolve().parents[1]
BASE='''<pf_dem>
<fluid dx=".1" rho0="1" sound_speed="10" viscosity=".01" sff=".001" sfs=".003"/>
<particle radius=".2" density="2" vx=".01" wz=".1"/>
<contacts stiffness="0" damping="0" friction="0"/>
<liquid mode="droplets"><droplet radius=".3" z=".25"/></liquid>
<domain margin=".1" max_samples="10000"/>
<run end_time="4e-5" dt_scale=".75" record_every_steps="1" dump_every_steps="0"/>
</pf_dem>'''
class FluidEquilibrationTests(unittest.TestCase):
 def test_xml_optional_runtime_setting(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'in.xml';p.write_text(BASE);a=parse_xml(p)
   self.assertEqual(a['run']['fluid_equilibration_time'],0)
   for value in ['0','2e-5']:
    p.write_text(BASE.replace('<run ','<run fluid_equilibration_time="'+value+'" '));b=parse_xml(p)
    self.assertEqual(b['run']['fluid_equilibration_time'],float(value))
    self.assertEqual(a['signature'],b['signature']);self.assertEqual(prepared_text(a),prepared_text(b))
   for value in ['-1','nan','inf']:
    p.write_text(BASE.replace('<run ','<run fluid_equilibration_time="'+value+'" '))
    with self.assertRaises(ValueError):parse_xml(p)
 def test_wrapper_forwards_seconds(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td);p=td/'in.xml';p.write_text(BASE.replace('<run ','<run fluid_equilibration_time="2e-5" '))
   r=subprocess.run([sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(p),'--output',str(td/'out'),'--prepare-only'],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr)
   cmd=json.loads((td/'out/command.json').read_text());i=cmd.index('--fluid-equilibration-time')
   self.assertEqual(float(cmd[i+1]),2e-5)
 def test_ke_option_is_runtime_only(self):
  with tempfile.TemporaryDirectory() as td:
   xml=Path(td)/'in.xml'
   xml.write_text(BASE)
   baseline=parse_xml(xml)
   xml.write_text(BASE.replace('<run ','<run ke_target_fraction=".01" ke_hold_time="2e-5" '))
   config=parse_xml(xml)
   self.assertEqual(config['signature'],baseline['signature'])
   self.assertEqual(prepared_text(config),prepared_text(baseline))
   self.assertEqual(config['run']['ke_target_fraction'],.01)
   self.assertEqual(config['run']['ke_hold_time'],2e-5)
   for value in ('0','1','nan'):
    xml.write_text(BASE.replace('<run ',f'<run ke_target_fraction="{value}" '))
    with self.assertRaises(ValueError):parse_xml(xml)
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for native kinetic-energy stop smoke')
 def test_ke_stop_after_release_with_final_outputs_and_time_cap(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td)
   xml=td/'in.xml'
   xml.write_text(BASE.replace('end_time="4e-5"','end_time=".4" fluid_equilibration_time=".05" density_check="false"')
                 .replace('viscosity=".01"','viscosity=".1"')
                 .replace('record_every_steps="1" dump_every_steps="0"',
                          'record_every_steps="99999" dump_every_steps="99999"'))
   def run(name,fraction,hold):
    out=td/name
    result=subprocess.run([sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(xml),
         '--executable',os.environ['PF_DEM_TEST_EXE'],'--output',str(out),'--end-time','.5',
         '--ke-target-fraction',fraction,'--ke-hold-time',hold],capture_output=True,text=True)
    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
    meta=json.loads((out/'dem_run.json').read_text())
    history=list(csv.DictReader((out/'dem_history.csv').read_text().splitlines()))
    self.assertEqual([row['fluid_equilibrating'] for row in history],['1','0','0'])
    self.assertEqual(float(history[1]['time']),.05)
    self.assertEqual(float(history[-1]['time']),meta['time'])
    self.assertTrue((out/f"liquid_{meta['step']:010}.vtp").is_file())
    self.assertTrue((out/'dem_checkpoint.dat').is_file())
    self.assertAlmostEqual(float(history[-1]['kinetic_energy'])/meta['peak_kinetic_energy'],
                           meta['final_ke_fraction'])
    return meta
   crossing=run('crossing','.8','0')
   held=run('held','.8','.02')
   for meta in (crossing,held):
    self.assertEqual(meta['termination_reason'],'kinetic_energy')
    self.assertFalse(meta['completed_time'])
    self.assertLessEqual(meta['final_ke_fraction'],.8)
    self.assertTrue(.05<meta['time']<.4)
   self.assertGreaterEqual(held['time']-crossing['time'],.02)
   capped=run('capped','.01','1')
   self.assertEqual(capped['termination_reason'],'time_limit')
   self.assertTrue(capped['completed_time'])
   self.assertEqual(capped['time'],.4)
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for native event/restart smoke')
 def test_native_release_and_resume(self):
  exe=os.environ['PF_DEM_TEST_EXE']
  with tempfile.TemporaryDirectory() as td:
   td=Path(td);xml=td/'in.xml';xml.write_text(BASE)
   config=parse_xml(xml);prepared=td/'prepared.dat';prepared.write_text(prepared_text(config))
   def run(name,flags):
    cwd=td/name;cwd.mkdir()
    r=subprocess.run([exe,'--input',str(prepared)]+flags,cwd=cwd,capture_output=True,text=True)
    return cwd,r
   out,r=run('whole',['--fluid-equilibration-time','2e-5'])
   self.assertEqual(r.returncode,0,r.stderr)
   rows=list(csv.DictReader((out/'dem_history.csv').read_text().splitlines()))
   self.assertEqual([int(x['fluid_equilibrating']) for x in rows],[1,0,0])
   self.assertEqual(float(rows[1]['time']),2e-5)
   bodies=list(csv.DictReader((out/'dem_particles.csv').read_text().splitlines()))
   self.assertEqual([float(bodies[1][k]) for k in ['x','y','z']],[0,0,0])
   self.assertEqual(float(bodies[1]['vx']),.01)
   self.assertGreater(float(bodies[-1]['x']),0)
   part,r=run('part',['--fluid-equilibration-time','2e-5','--end-time','1e-5'])
   self.assertEqual(r.returncode,0,r.stderr)
   resumed,r=run('resumed',['--resume',str(part/'dem_checkpoint.dat')])
   self.assertEqual(r.returncode,0,r.stderr)
   history=list(csv.DictReader((resumed/'dem_history.csv').read_text().splitlines()))
   self.assertEqual([int(x['fluid_equilibrating']) for x in history],[1,0,0])
   self.assertEqual(float(history[1]['time']),2e-5)
   _,r=run('conflict',['--resume',str(part/'dem_checkpoint.dat'),'--fluid-equilibration-time','3e-5'])
   self.assertEqual(r.returncode,1);self.assertIn('differs from checkpoint',r.stderr)
   # A completed hold must stay completed on a subsequent invocation.
   after,r=run('after',['--resume',str(out/'dem_checkpoint.dat'),'--end-time','5e-5','--fluid-equilibration-time','2e-5'])
   self.assertEqual(r.returncode,0,r.stderr)
   history=list(csv.DictReader((after/'dem_history.csv').read_text().splitlines()))
   self.assertTrue(all(x['fluid_equilibrating']=='0' for x in history))
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for native dose/restart smoke')
 def test_dose_then_continue_hold_and_add_again(self):
  exe=os.environ['PF_DEM_TEST_EXE']
  with tempfile.TemporaryDirectory() as td:
   td=Path(td)
   xml=td/'in.xml'
   xml.write_text(BASE.replace('<run ','<run fluid_equilibration_time="1e-5" '))
   def run(name,end,previous=None,add=False,expected=0):
    out=td/name
    cmd=[sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(xml),
         '--executable',exe,'--output',str(out),'--end-time',str(end)]
    if previous:
     cmd+=['--resume',str(previous/'dem_checkpoint.dat')]
    if add:
     cmd+=['--add-liquid-volume','.0031','--post-addition-equilibration-time','2e-5']
    result=subprocess.run(cmd,capture_output=True,text=True)
    self.assertEqual(result.returncode,expected,result.stdout+result.stderr)
    return out
   def history(out):
    return list(csv.DictReader((out/'dem_history.csv').read_text().splitlines()))
   seed=run('seed',2e-5)
   mid=run('mid',3e-5,seed,add=True)
   run('double-during-hold',5e-5,mid,add=True,expected=1)
   relaxed=run('relaxed',5e-5,mid)
   second=run('second',8e-5,relaxed,add=True)
   for before,after in ((seed,mid),(relaxed,second)):
    prior=history(before)[-1]
    first=history(after)[0]
    event=json.loads((after/'liquid_addition.json').read_text())
    self.assertEqual(first['step'],prior['step'])
    self.assertEqual(first['fluid_equilibrating'],'1')
    self.assertAlmostEqual(float(first['liquid_volume'])-float(prior['liquid_volume']),.0031,places=14)
    self.assertAlmostEqual(event['added_volume_m3'],.0031,places=14)
    self.assertGreater(event['added_particles'],0)
   self.assertEqual(history(relaxed)[0]['fluid_equilibrating'],'1')
   self.assertEqual(history(relaxed)[-1]['fluid_equilibrating'],'0')
   self.assertEqual(history(second)[-1]['fluid_equilibrating'],'0')
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for native uniform dose smoke')
 def test_uniform_doses_on_selected_surfaces(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td)
   xml=td/'in.xml'
   xml.write_text(BASE.replace('<run ','<run fluid_equilibration_time="2e-5" ')
                  .replace('dump_every_steps="0"','dump_every_steps="1"'))
   def run(name,end,flags):
    out=td/name
    cmd=[sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(xml),
         '--executable',os.environ['PF_DEM_TEST_EXE'],'--output',str(out),
         '--end-time',end]+flags
    result=subprocess.run(cmd,capture_output=True,text=True)
    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
    return out
   def points(path):
    values=map(float,ET.parse(path).find('.//Points/DataArray').text.split())
    return list(zip(values,values,values))
   seed=run('seed','2e-5',[])
   old=points(sorted(seed.glob('liquid_*.vtp'))[-1])
   for surface in ['liquid','solid']:
    with self.subTest(surface=surface):
     mode=f'uniform-{surface}-vapor'
     out=run(surface,'3e-5',['--resume',str(seed/'dem_checkpoint.dat'),
         '--add-liquid-volume','.006','--liquid-addition-mode',mode])
     event=json.loads((out/'liquid_addition.json').read_text())
     self.assertEqual(event['addition_mode'],mode)
     self.assertAlmostEqual(event['added_volume_m3'],.006,places=14)
     combined=points(sorted(out.glob('liquid_*.vtp'))[0])
     self.assertEqual(combined[:len(old)],old)
     added=combined[len(old):]
     self.assertEqual(len(added),event['added_particles'])
     self.assertGreater(event['eligible_sites'],len(added))
     for x in added:
      wall_gap=math.dist(x,(0,0,0))-.2
      liquid_gap=min(math.dist(x,y) for y in old)
      self.assertGreaterEqual(wall_gap,.05-1e-14)
      if surface=='liquid':
       self.assertTrue(.085-1e-14 <= liquid_gap <= .145+1e-14)
      else:
       self.assertGreater(liquid_gap,.145)
       self.assertLessEqual(wall_gap,.15+1e-14)
 @unittest.skipUnless(os.environ.get('PF_DEM_TEST_EXE'),'Set PF_DEM_TEST_EXE for native Kelvin dose smoke')
 def test_multilayer_doses_and_kelvin_weights(self):
  with tempfile.TemporaryDirectory() as td:
   td=Path(td)
   xml=td/'in.xml'
   xml.write_text(BASE.replace('dump_every_steps="0"','dump_every_steps="1"'))
   def run(name,end,flags,expected=0):
    out=td/name
    cmd=[sys.executable,str(ROOT/'tools/run_pf_dem.py'),str(xml),
         '--executable',os.environ['PF_DEM_TEST_EXE'],'--output',str(out),
         '--end-time',end]+flags
    result=subprocess.run(cmd,capture_output=True,text=True)
    self.assertEqual(result.returncode,expected,result.stdout+result.stderr)
    return out,result.stderr
   seed,_=run('seed','2e-5',[])
   checkpoint=seed/'dem_checkpoint.dat'
   original=checkpoint.read_bytes()
   dose=['--resume',str(checkpoint),'--liquid-addition-mode','curvature-based']
   def insert(name,rng,nu='1.5',length='.02',expected=0,volume='.012',end='3e-5'):
    return run(name,end,dose+['--add-liquid-volume',volume,'--saturation-ratio',nu,
                             '--kelvin-length',length,'--liquid-addition-seed',rng],expected)
   first,_=insert('first','51')
   repeat,_=insert('repeat','51')
   different,_=insert('different','89')
   for out in (first,repeat,different):
    event=json.loads((out/'liquid_addition.json').read_text())
    self.assertEqual(event['layers'],1)
    self.assertAlmostEqual(event['added_volume_m3'],.012,places=14)
    self.assertGreater(event['eligible_liquid_vapor_sites'],0)
    self.assertGreater(event['eligible_solid_vapor_sites'],0)
    self.assertEqual(event['added_solid_vapor_particles']+event['added_liquid_vapor_particles'],
                     event['added_particles'])
   def snapshot(out):
    return sorted(out.glob('liquid_*.vtp'))[0].read_bytes()
   self.assertEqual(snapshot(first),snapshot(repeat))
   self.assertNotEqual(snapshot(first),snapshot(different))
   kelvin_layers,_=insert('kelvin-layers','42',volume='.55',end='2.00001e-5')
   liquid_vapor_layers,_=run('liquid-vapor-layers','2.00001e-5',[
       '--resume',str(checkpoint),'--liquid-addition-mode','liquid-vapor',
       '--add-liquid-volume','.55','--post-addition-equilibration-time','1e-5'])
   for multilayer,min_layers in ((kelvin_layers,3),(liquid_vapor_layers,2)):
    with self.subTest(mode=multilayer.name):
     event=json.loads((multilayer/'liquid_addition.json').read_text())
     self.assertGreaterEqual(event['layers'],min_layers)
     self.assertAlmostEqual(event['added_volume_m3'],.55,places=14)
     history=list(csv.DictReader((multilayer/'dem_history.csv').read_text().splitlines()))
     self.assertEqual(float(history[0]['time']),2e-5)
     self.assertAlmostEqual(float(history[0]['liquid_volume'])-event['volume_before_m3'],.55,places=14)
     points=ET.fromstring(snapshot(multilayer)).find('.//Points/DataArray').text.split()
     new=points[-3*event['added_particles']:]
     self.assertEqual(len(set(zip(new[::3],new[1::3],new[2::3]))),event['added_particles'])
     if multilayer==liquid_vapor_layers:
      self.assertEqual(history[0]['fluid_equilibrating'],'1')
   blocked_solid,_=insert('solid-cutoff','51',nu='1.5',length='.2')
   cutoff=json.loads((blocked_solid/'liquid_addition.json').read_text())
   self.assertGreater(cutoff['resolved_solid_vapor_sites'],0)
   self.assertEqual(cutoff['eligible_solid_vapor_sites'],0)
   self.assertGreater(cutoff['added_liquid_vapor_particles'],0)
   refused,error=insert('refused','51',nu='1',length='0',expected=1)
   self.assertIn('non-positive Kelvin weights',error)
   self.assertFalse((refused/'dem_checkpoint.dat').exists())
   self.assertEqual(checkpoint.read_bytes(),original)
if __name__=='__main__':unittest.main()
