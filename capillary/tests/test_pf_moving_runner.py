import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

P = Path(__file__).resolve().parents[1] / 'tools' / 'run_pf_moving_smoke.py'
spec = importlib.util.spec_from_file_location('run_pf_moving_smoke', P)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RunnerTest(unittest.TestCase):
    def test_short_output_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            (out / 'pf_checkpoint.dat').write_text('test')
            with (out / 'pf_history.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['time','linear_error','angular_error','max_density_error'])
                writer.writeheader()
                writer.writerows([{'time': '0', 'linear_error': '0', 'angular_error': '0', 'max_density_error': '0'},
                                  {'time': '.002', 'linear_error': '1e-12', 'angular_error': '2e-11', 'max_density_error': '.01'}])
            with (out / 'pf_rigid.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['step','cx','cy'])
                writer.writeheader()
                writer.writerows([{'step':'0','cx':'0','cy':'0'},
                                  {'step':'2','cx':'0.0001','cy':'0'}])
            report=module.check_outputs(out,.002,.05)
            self.assertGreater(report['final_x'],report['initial_x'])
            self.assertEqual(report['frames'],2)
            with self.assertRaises(RuntimeError):
                module.check_outputs(out,.004,.05)


if __name__ == '__main__':
    unittest.main()
