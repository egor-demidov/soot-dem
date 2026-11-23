"""Original LB XML structure mapped to the actual soot-dem/SPH driver."""
from pathlib import Path
import math,shutil,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from run_pf_dem import parse_xml,prepared_text,SOOT_DEFAULTS
ROOT=Path(__file__).resolve().parents[1]
class LegacyInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        for f in (ROOT/'inputs').glob('*.xml'):shutil.copyfile(f,self.p/f.name)
        aggregate=self.p/'aggregate_bank/aggregate_1';aggregate.mkdir(parents=True)
        (aggregate/'final_aggregate.vtk').write_text('# vtk DataFile Version 3.0\ntest\nASCII\nDATASET POLYDATA\nPOINTS 3 float\n-2 0 0 0 0 0 2 0 0\n')
    def tearDown(self):self.tmp.cleanup()
    def parse(self,name='capillary_levelset.xml'):return parse_xml(self.p/name)
    def test_original_geometry_units_and_dem_coefficients(self):
        a=self.parse();b=self.parse('capillary_dem.xml');m=a['model']
        self.assertEqual(m['bodies'],b['model']['bodies'])
        self.assertAlmostEqual(m['bodies'][0]['center'][0]/1e-9,-28)
        self.assertAlmostEqual(m['bodies'][0]['radius']/1e-9,14)
        self.assertAlmostEqual(m['fluid'][0]/1e-9,1.75)
        self.assertEqual(m['fluid'][1:4],[1258,2000,.001])
        self.assertAlmostEqual(m['rm']/1e-9,14)
        self.assertEqual((m['alpha'],m['threshold']),(.5,.5))
        values=dict(zip((k for k,_ in SOOT_DEFAULTS),m['soot_dem']['values']))
        self.assertEqual(values['k_n'],10000);self.assertEqual(values['gamma_n'],5e-9)
        self.assertEqual(values['A'],1e-19);self.assertEqual(values['max_dt'],5e-13)
        self.assertEqual(m['soot_dem']['neck_mode'],1)
        d=b['model']['drops'][0]
        for x,y in zip(d['center'],[-140.5052e-9,90.9034e-9,143.8088e-9]):self.assertAlmostEqual(x/y,1)
        self.assertAlmostEqual(d['radius']/1e-9,60)
        self.assertTrue(prepared_text(a).startswith('PFDEM_INPUT 2\n'))
    def test_targets_are_metadata_not_unverified_calibration(self):
        a=self.parse();path=self.p/'capillary_levelset.xml'
        path.write_text(path.read_text().replace('>60</let>','>120</let>'))
        b=self.parse();self.assertEqual(a['signature'],b['signature'])
        self.assertEqual(b['input_compatibility']['target_physics']['liquid_contact_angle_deg'],120)
    def test_fluid_equilibration_in_sph_or_let(self):
        a=self.parse();path=self.p/'capillary_levelset.xml';base=path.read_text()
        self.assertEqual(a['run']['fluid_equilibration_time'],0)
        for content in [base.replace('<sph ','<sph fluid_equilibration_time="1e-9" '),
                        base.replace('<initialization>','<let id="fluid_equilibration_time" type="real">1e-9</let><initialization>')]:
            path.write_text(content);b=self.parse()
            self.assertEqual(b['run']['fluid_equilibration_time'],1e-9)
            self.assertEqual(a['signature'],b['signature'])
        path.write_text(base.replace('<sph ','<sph fluid_equilibration_time="-1" '))
        with self.assertRaisesRegex(ValueError,'fluid_equilibration_time'):self.parse()
        path.write_text(base.replace('<sph ','<sph fluid_equilibration_time="1" ').replace('<initialization>','<let id="fluid_equilibration_time" type="real">1e-9</let><initialization>'))
        with self.assertRaisesRegex(ValueError,'Specify fluid_equilibration_time once'):self.parse()
    def test_missing_explicit_sph_fails(self):
        import xml.etree.ElementTree as ET
        path=self.p/'capillary_levelset.xml';tree=ET.parse(path);tree.getroot().remove(tree.getroot().find('sph'));tree.write(path)
        with self.assertRaisesRegex(ValueError,'explicit <sph'):self.parse()
    def test_neck_strengths_convert_to_native_force_thresholds(self):
        path=self.p/'capillary_levelset.xml'
        base=path.read_text()
        declarations={
            'neck_radius':'8e-9','normal_strength':'16e9','shear_strength':'12e9'}
        lets={key:f'<let id="{key}" type="real">{value}</let>'
              for key,value in declarations.items()}
        full_input=base.replace('<initialization>',''.join(lets.values())+'<initialization>')
        path.write_text(full_input)
        result=self.parse()
        coefficients=dict(zip((key for key,_ in SOOT_DEFAULTS),
                              result['model']['soot_dem']['values']))
        self.assertTrue(math.isclose(coefficients['normal_yield'],math.pi*(8e-9)**2*16e9))
        self.assertTrue(math.isclose(coefficients['shear_yield'],math.pi*(8e-9)**2*12e9))
        self.assertTrue(prepared_text(result).startswith('PFDEM_INPUT 2\n'))
        path.write_text(full_input.replace(lets['shear_strength'],''))
        with self.assertRaisesRegex(ValueError,'Specify neck_radius'):self.parse()
        path.write_text(full_input)
        necks=self.p/'capillary_necks.xml'
        necks.write_text(necks.read_text().replace('id="normal_yield" type="real">0',
                                                  'id="normal_yield" type="real">1'))
        with self.assertRaisesRegex(ValueError,'nonzero normal_yield'):self.parse()
    def test_includes_cycles_and_duplicates_rejected(self):
        path=self.p/'capillary_contacts.xml'
        path.write_text(path.read_text().replace('</include_file>','<include>capillary_contacts.xml</include></include_file>'))
        with self.assertRaisesRegex(ValueError,'Cyclic'):self.parse()
    def test_legacy_aliases_preserve_model_and_cadence(self):
        a=self.parse();path=self.p/'capillary_levelset.xml';s=path.read_text()
        s=s.replace('particles_per_primary_diameter="16"','').replace('<include>capillary_contacts.xml</include>','<let id="lb_cells_per_particle_diameter" type="integer">16</let><include>capillary_contacts.xml</include>')
        for x in ('surface','fluid','particle','log'):s=s.replace(x+'_every_sph_steps',x+'_every_lb_steps')
        s=s.replace('dem_every_sph_steps','dem_every_dem_steps');path.write_text(s)
        b=self.parse();self.assertEqual(a['signature'],b['signature']);self.assertEqual(a['run'],b['run'])
    def test_anisotropic_margins_and_unknown_settings(self):
        a=self.parse('capillary_dem.xml');path=self.p/'capillary_dem.xml'
        path.write_text(path.read_text().replace('x="30" y="30" z="30"','x="40" y="30" z="20"'))
        b=self.parse('capillary_dem.xml')
        self.assertAlmostEqual((a['model']['lower'][0]-b['model']['lower'][0])/1e-9,10)
        self.assertAlmostEqual((a['model']['upper'][2]-b['model']['upper'][2])/1e-9,10)
        path.write_text(path.read_text().replace('<initialization>','<let id="typo" type="real">1</let><initialization>'))
        with self.assertRaisesRegex(ValueError,'Unsupported let'):self.parse('capillary_dem.xml')
if __name__=='__main__':unittest.main()
