import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'tools'))
from run_pf_wetting_checks import commands

class WettingRunnerTests(unittest.TestCase):
    def test_two_matched_step_runs_preserve_geometry_and_strengths(self):
        plans=commands('pf_capillary_validation')
        self.assertEqual(len(plans),2)
        def opts(cmd):return dict(zip(cmd[1::2],cmd[2::2]))
        a,b=[opts(cmd) for _,cmd in plans]
        for key in ['--case','--dx','--radius','--sff','--sfs','--viscosity','--end-time']:
            self.assertEqual(a[key],b[key])
        self.assertEqual(a['--case'],'sessile')
        self.assertEqual(float(a['--sfs'])/float(a['--sff']),3.5)
        acoustic=[float(v['--dt-scale'])/float(v['--sound-speed']) for v in (a,b)]
        self.assertAlmostEqual(*acoustic,places=15)
        self.assertGreater(float(a['--end-time']),5.)
        for v in (a,b):
            self.assertGreater(float(v['--snapshot-interval']),0)
            self.assertEqual(v['--dump-every'],'0')

if __name__=='__main__':unittest.main()
