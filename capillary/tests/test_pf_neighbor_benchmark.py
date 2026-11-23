"""Check that A/B neighbor benchmark comparisons do not silently reorder particles."""
import csv
import importlib.util
import math
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[1] / 'tools' / 'benchmark_pf_neighbors.py'
spec = importlib.util.spec_from_file_location('benchmark_pf_neighbors', MODULE)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def write_state(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['x', 'y', 'z', 'rho', 'vx', 'vy', 'vz'])
        writer.writerows(data)


def test_state_comparison(tmp_path):
    original = [[0., .1, -.1, 1., 0., 0., 0.], [.5, .2, -.2, 1.01, .1, .2, .3]]
    changed = [[0., .1, -.1, 1., 0., 0., 0.], [.500001, .2, -.2, 1.01, .1, .2, .3]]
    write_state(tmp_path / 'native' / 'pf_final_liquid.csv', original)
    write_state(tmp_path / 'verlet' / 'pf_final_liquid.csv', changed)
    result = benchmark.compare_state(tmp_path)
    assert result['particles'] == 2
    assert math.isclose(result['max_abs_difference']['x'], 1e-6, rel_tol=1e-9)
    assert result['max_abs_difference']['rho'] == 0.


def test_comparison_rejects_missing_particle(tmp_path):
    write_state(tmp_path / 'native' / 'pf_final_liquid.csv', [[0] * 7])
    write_state(tmp_path / 'verlet' / 'pf_final_liquid.csv', [[0] * 7, [1] * 7])
    with pytest.raises(RuntimeError, match='counts mismatch'):
        benchmark.compare_state(tmp_path)
