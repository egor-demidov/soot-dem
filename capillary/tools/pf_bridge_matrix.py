#!/usr/bin/env python3
"""Prepare, run, and fit a 2×2 PF-SPH fixed-sphere bridge matrix."""
import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from young_laplace import bridge_reference


@dataclass(frozen=True)
class Case:
    radius_nm: float
    gap_nm: float

    @property
    def name(self):
        return f"r{self.radius_nm:g}-g{self.gap_nm:g}"


CASES = [Case(radius, gap) for radius in (21., 31.5) for gap in (7., 10.5)]
HOLDOUT = Case(31.5, 7.)


def prepare(root, sfs, end_time, threads):
    root.mkdir(parents=True, exist_ok=False)
    for case in CASES:
        center = (case.radius_nm + case.gap_nm / 2) * 1e-9
        radius = case.radius_nm * 1e-9
        seed = radius * 5 / 6
        (root / f"{case.name}.xml").write_text(f"""<pf_dem>
  <fluid dx="1.75e-9" rho0="1258" sound_speed="1000"
         viscosity="0.001" sff="0.002" sfs="{sfs}"/>
  <particle x="{-center:.12g}" y="0" z="0" radius="{radius:.12g}" density="1700" fixed="true"/>
  <particle x="{center:.12g}" y="0" z="0" radius="{radius:.12g}" density="1700" fixed="true"/>
  <liquid mode="droplets">
    <droplet x="0" y="0" z="0" radius="{seed:.12g}"/>
  </liquid>
  <run end_time="{end_time:.12g}" dt_scale="0.75" threads="{threads}"
       record_every_steps="100" dump_every_steps="1000"
       checkpoint_every_steps="1000"/>
</pf_dem>
""")
        print(root / f"{case.name}.xml")


def run_matrix(root, executable):
    runner = Path(__file__).with_name("run_pf_dem.py")
    failures = []
    for case in CASES:
        source = root / f"{case.name}.xml"
        if not source.is_file():
            raise ValueError(f"Missing {source}. Run prepare first.")
        out = root / case.name
        if out.exists():
            status = out / "dem_run.json"
            if status.is_file() and json.loads(status.read_text()).get("completed_time"):
                print(f"Using completed case {out}", flush=True)
                continue
            print(f"Existing incomplete case {out}. Resume it or move the directory aside.", flush=True)
            failures.append(case.name)
            continue
        result = subprocess.run([
            sys.executable, str(runner), str(source),
            "--executable", executable, "--output", str(out),
        ])
        if result.returncode:
            failures.append(case.name)
    if failures:
        raise ValueError("Incomplete cases: " + ", ".join(failures))


def read_case(root, case, gamma):
    folder = root / case.name
    model = json.loads((folder / "resolved_input.json").read_text())["model"]
    status = json.loads((folder / "dem_run.json").read_text())
    if not status.get("completed_time"):
        raise ValueError(f"{case.name} has not completed")
    if model.get("wall") or len(model["bodies"]) != 2:
        raise ValueError(f"{case.name}: expected two spheres and no wall")
    bodies = model["bodies"]
    if not all(b["fixed"] and math.isclose(b["radius"], case.radius_nm * 1e-9) for b in bodies):
        raise ValueError(f"{case.name}: incorrect sphere radius or fixity")
    offset = np.subtract(bodies[1]["center"], bodies[0]["center"])
    distance = np.linalg.norm(offset)
    gap = float(distance - 2 * bodies[0]["radius"])
    if not math.isclose(gap, case.gap_nm * 1e-9, rel_tol=1e-6):
        raise ValueError(f"{case.name}: incorrect surface gap")
    axis = offset / distance
    history = np.genfromtxt(folder / "dem_history.csv", delimiter=",", names=True, ndmin=1)
    particles = np.genfromtxt(folder / "dem_particles.csv", delimiter=",", names=True, ndmin=1)
    t = history["time"]
    late = history[t >= t[-1] - .25 * (t[-1] - t[0])]
    if len(late) < 5 or np.any(np.diff(late["time"]) <= 0):
        raise ValueError(f"{case.name}: need at least five distinct late history rows")
    axial = []
    for i in (0, 1):
        p = particles[(particles["particle"] == i) & np.isin(particles["step"], late["step"])]
        if not np.array_equal(p["step"], late["step"]):
            raise ValueError(f"{case.name}: particle loads do not match the history steps")
        axial.append(np.column_stack([p[f"fluid_f{k}"] for k in "xyz"]) @ axis)
    a, b = axial
    force = (a - b) / 2
    radius = bodies[0]["radius"]
    scale = gamma * radius
    mean = float(force.mean())
    volume = float(late["liquid_volume"].mean())
    mach = float(late["max_speed"].max() / model["fluid"][2])
    force_range = float(np.ptp(force))
    imbalance = float(np.max(abs(a+b)))
    if not (np.isfinite(force).all() and math.isfinite(mach) and math.isfinite(volume) and volume > 0):
        raise ValueError(f"{case.name}: invalid late force, volume or Mach")
    force_tolerance = .1 * max(abs(mean), .1 * scale)
    checks = []
    if force_range > force_tolerance:
        checks.append("late force variation")
    if imbalance > force_tolerance:
        checks.append("axial force imbalance")
    if mach > .02:
        checks.append("late Mach above 0.02")
    if float(np.ptp(late["liquid_volume"])) > .01 * volume:
        checks.append("late liquid volume variation")
    return dict(
        case=case.name, radius_m=radius, gap_m=gap, volume_m3=volume,
        sph_force_N=mean, late_force_range_N=force_range,
        max_axial_imbalance_N=imbalance, late_mach=mach,
        late_rows=len(late), screening_failures=checks, fluid=model["fluid"],
    )


def analyze(root, gamma, limits):
    cases = [read_case(root, case, gamma) for case in CASES]
    fluid = cases[0]["fluid"]
    if any(r["fluid"] != fluid for r in cases):
        raise ValueError("Use identical fluid parameters and sfs for all four cases")
    if any(r["screening_failures"] for r in cases):
        bad = {r["case"]: r["screening_failures"] for r in cases if r["screening_failures"]}
        raise ValueError(f"Relaxation screen failed: {bad}")

    @lru_cache(maxsize=None)
    def reference_force(index, angle):
        row = cases[index]
        try:
            return bridge_reference(
                row["volume_m3"], row["gap_m"], angle,
                radius=row["radius_m"], sigma=gamma,
            )["attractive_force"]
        except RuntimeError:
            return np.nan

    def discrepancy(i, angle):
        row = cases[i]
        return (reference_force(i, float(angle)) - row["sph_force_N"]) / (gamma * row["radius_m"])

    n = math.ceil((limits[1] - limits[0]) / 10)
    angles = np.linspace(*limits, n + 1)
    for i, row in enumerate(cases):
        errors = [discrepancy(i, angle) for angle in angles]
        roots = [float(t) for t, error in zip(angles, errors) if abs(error) < 1e-10]
        for j in range(n):
            if errors[j] * errors[j + 1] < 0:
                try:
                    root = float(brentq(lambda t: discrepancy(i, t),
                                        angles[j], angles[j + 1], xtol=1e-5))
                    if abs(discrepancy(i, root)) < 1e-5:
                        roots.append(root)
                except ValueError:
                    # A reference branch disappeared within the bracket.
                    continue
        row["force_equivalent_angles_deg"] = sorted(roots)
        if not roots:
            print(f"{row['case']}: no force crossing found in the angle scan", file=sys.stderr)

    def fit(indices):
        def loss(angle):
            errors = np.array([discrepancy(i, angle) for i in indices])
            return float(np.mean(errors**2)) if np.isfinite(errors).all() else np.inf

        values = np.array([loss(angle) for angle in angles])
        if not np.isfinite(values).any():
            raise ValueError("No common Young–Laplace branch found in the angle scan")
        j = int(values.argmin())
        best = float(angles[j])
        opt = minimize_scalar(loss, method="bounded", options={"xatol": 1e-4},
                              bounds=(angles[max(0, j - 1)], angles[min(n, j + 1)]))
        if opt.success and opt.fun < values[j]:
            best = float(opt.x)
        return best, math.sqrt(loss(best))

    held = CASES.index(HOLDOUT)
    common, rms = fit(tuple(range(len(cases))))
    trained, _ = fit(tuple(i for i in range(len(cases)) if i != held))
    predicted = reference_force(held, trained)
    for i, row in enumerate(cases):
        prediction = reference_force(i, common)
        row["common_angle_force_N"] = float(prediction)
        row["common_angle_force_error_N"] = float(prediction - row["sph_force_N"])
    return dict(
        sfs=fluid[5], gamma_N_per_m=gamma, angle_scan_deg=limits,
        common_angle_deg=common, rms_force_error_over_gamma_R=rms,
        common_angle_at_scan_limit=bool(min(abs(common - t) for t in limits) < 1e-3),
        holdout=dict(case=HOLDOUT.name, fitted_on_other_three_deg=trained,
                     sph_force_N=cases[held]["sph_force_N"],
                     predicted_force_N=float(predicted) if np.isfinite(predicted) else None),
        cases=cases,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("matrix", type=Path)
    p.add_argument("--sfs", type=float, default=.008)
    p.add_argument("--end-time", type=float, default=2e-9)
    p.add_argument("--threads", type=int, default=4)
    r = commands.add_parser("run")
    r.add_argument("matrix", type=Path)
    r.add_argument("--executable", default="pf_soot_dem")
    a = commands.add_parser("analyze")
    a.add_argument("matrix", type=Path)
    a.add_argument("--gamma", type=float, default=.08781)
    a.add_argument("--theta-range", type=float, nargs=2, default=(15., 165.))
    args = parser.parse_args()
    try:
        if args.mode == "prepare":
            if not (math.isfinite(args.sfs) and args.sfs >= 0 and
                    math.isfinite(args.end_time) and args.end_time > 0 and args.threads > 0):
                raise ValueError("Require finite sfs ≥ 0, positive time and positive thread count")
            prepare(args.matrix, args.sfs, args.end_time, args.threads)
        elif args.mode == "run":
            run_matrix(args.matrix, args.executable)
        else:
            if not (math.isfinite(args.gamma) and args.gamma > 0
                    and 0 < args.theta_range[0] < args.theta_range[1] < 180):
                raise ValueError("Require positive gamma and increasing angle limits in (0, 180)")
            report = analyze(args.matrix, args.gamma, args.theta_range)
            print(json.dumps(report, indent=2, allow_nan=False))
    except (ValueError, OSError, KeyError, IndexError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
