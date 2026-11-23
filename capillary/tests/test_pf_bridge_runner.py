"""Dry-run command and validation logic of the native two-sphere bridge runner."""
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('bridge_runner', ROOT / 'tools' / 'run_pf_bridge_smoke.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_bridge_dry_run():
    cmd = [sys.executable, str(ROOT / 'tools' / 'run_pf_bridge_smoke.py'),
           '--dry-run', '--executable', 'fake.exe', '--fixed-spheres', '0',
           '--gap', '0.25', '--end-time', '0.001']
    out = subprocess.run(cmd, text=True, capture_output=True)
    assert out.returncode == 0, out.stderr
    assert '--case bridge' in out.stdout
    assert '--fixed-spheres 0' in out.stdout
    assert '--bridge-gap 0.25' in out.stdout


def test_bridge_reject_invalid():
    for args in (['--gap', '-0.1'], ['--dx', '0'], ['--fixed-spheres', '8'],
                 ['--contact-stiffness', '-1']):
        run = subprocess.run([sys.executable, str(ROOT / 'tools' / 'run_pf_bridge_smoke.py'),
                              '--dry-run', *args], capture_output=True, text=True)
        assert run.returncode != 0


def test_assess_reject_missing(tmp_path):
    from types import SimpleNamespace
    try:
        module.assess(tmp_path, SimpleNamespace())
    except RuntimeError as exc:
        assert 'Missing output' in str(exc)
    else:
        raise AssertionError('Missing bridge outputs must be rejected')
