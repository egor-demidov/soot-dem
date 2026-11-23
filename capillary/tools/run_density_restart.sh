#!/usr/bin/env bash
# Short, read-only restart investigation. Never overwrites the source run.
set -u
if [ "$#" -lt 3 ]; then
  echo "Usage: bash $0 INPUT.xml CHECKPOINT.dat PF_SOOT_DEM [NEW_OUTPUT_DIR]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
python3 "$script_dir/run_pf_dem.py" "$1" --resume "$2" --executable "$3" \
  --output "${4:-results-density-debug}" --end-time 3.4e-10 --steps 5000
run_status=$?
python3 - "${4:-results-density-debug}" <<'PY'
import csv,sys
from pathlib import Path
from collections import deque
path=Path(sys.argv[1])/'density_extrema.csv'
if not path.is_file():
    print('No density diagnostics found; inspect run.log and ensure the executable was rebuilt.')
    sys.exit(1)
with path.open() as f:rows=deque(csv.DictReader(f),maxlen=12)
for r in rows:
    print(f"t_ns={float(r['time'])*1e9:.8g} {r['extreme']} sample={r['sample']} "
          f"rho/rho0={float(r['rho_over_rho0']):.8g} "
          f"P_Pa={float(r['pressure_pa']):.8g} "
          f"div_s^-1={float(r['continuity_divergence']):.8g} "
          f"body={r['nearest_body']} gap_nm={float(r['surface_gap'])*1e9:.8g} "
          f"neighbors_liquid={r['liquid_neighbors']} boundary={r['boundary_neighbors']}")
PY
exit "$run_status"
