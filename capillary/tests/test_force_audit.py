import csv
import pathlib
import sys
import tempfile
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'tools'))
from force_audit import force_audit

class ForceAuditTests(unittest.TestCase):
    def test_external_force_accounting_retains_raw_defect(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp)
            keys=['p'+a for a in 'xyz']+[p+a for p in ('support_impulse_','capillary_impulse_') for a in 'xyz']
            rows=[dict.fromkeys(keys,0.) for _ in range(2)]
            rows[1]['support_impulse_z']=-3.
            rows[1]['capillary_impulse_z']=-2.
            with (root/'history.csv').open('w') as f:
                w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
            result=force_audit(root)
            self.assertEqual(result['raw_max_supported_momentum_defect'],3.)
            self.assertEqual(result['max_residual_treating_CSF_as_external'],1.)
            self.assertEqual(result['final_supported_momentum_defect'],[0.,0.,-3.])
            self.assertIn('diagnostic only',result['warning'])
    def test_old_history_is_supported(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp);(root/'history.csv').write_text('time\n0\n')
            self.assertIsNone(force_audit(root))

if __name__=='__main__': unittest.main()
