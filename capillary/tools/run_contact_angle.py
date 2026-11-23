#!/usr/bin/env python3
"""Two focused CA trajectories. No builds, parameter sweeps, or VTP dumps."""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess
import time
import zipfile
from analyze import analyze
ROOT = pathlib.Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', default='capillary_validation')
    p.add_argument('--output', type=pathlib.Path, default=ROOT/'results-ca-long')
    p.add_argument('--angles', nargs='+', type=float, default=[60.,120.])
    p.add_argument('--end-time', type=float, default=40.)
    p.add_argument('--dx', type=float, default=.06)
    p.add_argument('--threads', type=int, default=6)
    p.add_argument('--snapshot-interval', type=float, default=1.)
    p.add_argument('--timeout', type=float, default=172800, help='Wall seconds per case')
    p.add_argument('--viscosity-update', choices=['acoustic','advection'], default='acoustic')
    p.add_argument('--audit-every', type=int, default=0, help='Momentum stages and gas outliers; 0 disables')
    p.add_argument('--smoke', action='store_true', help='Only two acoustic steps; NOT physical validation')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    if not (0<a.dx<=.3 and a.end_time>0 and a.threads>0 and a.snapshot_interval>0 and a.timeout>0 and a.audit_every>=0):
        p.error('Positive time, spacing, interval, timeout and thread count required; dx <= .3')
    if not 1<=len(a.angles)<=2 or len(set(a.angles))!=len(a.angles) or any(not 0<t<180 for t in a.angles):
        p.error('Choose one or two distinct angles between 0 and 180 degrees')
    exe = shutil.which(a.executable)
    if not exe and pathlib.Path(a.executable).is_file(): exe = str(pathlib.Path(a.executable).resolve())
    if not exe and not a.dry_run: p.error('Executable not found; rebuild first or specify --executable')
    exe = str(pathlib.Path(exe).resolve()) if exe else a.executable
    out = a.output.resolve()
    entries=[]
    for theta in a.angles:
        name=f'sessile_dx{.3 if a.smoke else a.dx:g}_theta{theta:g}'
        cmd=[exe,'--case','sessile','--dx',str(.3 if a.smoke else a.dx),'--theta',str(theta),
             '--end-time',str(.01 if a.smoke else a.end_time),'--steps',str(2 if a.smoke else 10000000),
             '--threads',str(a.threads),'--record-every',str(1 if a.smoke else 100),
             '--dump-every','0','--snapshot-interval',str(a.snapshot_interval),
             '--viscosity-update',a.viscosity_update]
        if a.audit_every: cmd += ['--audit-every',str(a.audit_every)]
        entries.append({'name':name,'command':cmd})
    if a.dry_run:
        print(json.dumps(entries,indent=2)); return 0
    if '--snapshot-interval' not in subprocess.check_output([exe,'--help'],text=True):
        p.error('Executable is old: rebuild capillary_validation with updated sources')
    if a.audit_every and '--audit-every' not in subprocess.check_output([exe,'--help'],text=True):
        p.error('Executable is old: rebuild for --audit-every')
    if '--viscosity-update' not in subprocess.check_output([exe,'--help'],text=True):
        p.error('Executable is old: rebuild for the viscosity scheduling update')
    if out.exists(): p.error('Output exists; choose a fresh --output')
    out.mkdir(parents=True)
    manifest={'smoke_only':a.smoke,'cases':entries,'executable_sha256':hashlib.sha256(pathlib.Path(exe).read_bytes()).hexdigest(),
              'sources':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
                         for sub in ('src','tools') for f in (ROOT/sub).rglob('*') if f.suffix in ('.h','.cpp','.py')}}
    failed=False
    try:
        for entry in entries:
            dest=out/entry['name']; dest.mkdir(); start=time.monotonic()
            print('Running',entry['name'],'(no VTP output)',flush=True)
            try:
                with (dest/'run.log').open('w') as log:
                    proc=subprocess.run(entry['command'],cwd=dest,stdout=log,stderr=subprocess.STDOUT,timeout=a.timeout)
                entry['returncode']=proc.returncode
                if proc.returncode==0:
                    entry['analysis']=analyze(dest)
                    if not a.smoke and not entry['analysis']['execution_completed']: failed=True
                else: failed=True
            except Exception as e:
                entry['error']=str(e); failed=True
            entry['wall_seconds']=time.monotonic()-start
            (out/'summary.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    finally:
        (out/'summary.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
        archive=out.with_suffix('.zip')
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
            for f in out.rglob('*'):
                if f.is_file(): z.write(f,f.relative_to(out.parent))
        print('Send back:',archive,flush=True)
    return int(failed)

if __name__=='__main__':
    raise SystemExit(main())
