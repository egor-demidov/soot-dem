#!/usr/bin/env python3
"""Two longer, matched-acoustic-step wetting runs; no parameter sweep or free-droplet rerun."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import zipfile
from run_pf_validation import analyze


def commands(executable,end_time=15.,threads=6,snapshots=1.,dump=0):
    common=['--case','sessile','--dx','.06','--radius','.6','--sff','.001',
            '--sfs','.0035','--viscosity','.05','--end-time',str(end_time),
            '--threads',str(threads),'--snapshot-interval',str(snapshots),
            '--dump-every',str(dump),'--record-every','100']
    return [(name,[executable]+common+['--sound-speed',str(sound),'--dt-scale',str(scale)])
            for name,sound,scale in [('c10_dt_half',10.,.5),('c20_dt_full',20.,1.)]]


def inspect(folder,grid_dx):
    summary=analyze(folder)
    from pf_interface import analyze_folder
    summary['interface']=analyze_folder(folder,grid_dx)
    (folder/'wetting_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',default='pf_capillary_validation')
    p.add_argument('--output',type=Path,default=Path('results-pf-wetting-checks'))
    p.add_argument('--end-time',type=float,default=15.)
    p.add_argument('--threads',type=int,default=6)
    p.add_argument('--snapshot-interval',type=float,default=1.)
    p.add_argument('--grid-dx',type=float,default=.25)
    p.add_argument('--dump-every',type=int,default=0)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--analyze-only',action='store_true',help='Reanalyze existing output; run no executable')
    a=p.parse_args()
    if not (math.isfinite(a.end_time) and a.end_time>0 and a.threads>0 and
            math.isfinite(a.snapshot_interval) and a.snapshot_interval>0 and
            .05<=a.grid_dx<=1 and a.dump_every>=0):
        p.error('Invalid time, spacing, output interval or thread count')
    exe=shutil.which(a.executable)
    if not exe and Path(a.executable).is_file():exe=str(Path(a.executable).resolve())
    if not exe and not a.dry_run and not a.analyze_only:p.error('Executable not found; rebuild first')
    plans=commands(exe or a.executable,a.end_time,a.threads,a.snapshot_interval,a.dump_every)
    if a.dry_run:
        for name,cmd in plans:print(name+': '+shlex.join(cmd))
        return
    # Check offline dependencies BEFORE spending time on simulations.
    import numpy,scipy,skimage.measure  # noqa: F401
    from pf_interface import analyze_folder  # noqa: F401
    archive=a.output.resolve().with_suffix('.zip')
    if a.analyze_only:
        if not a.output.is_dir():p.error('Existing output directory required')
    else:
        if a.output.exists() or archive.exists():p.error('Choose a new output directory/archive name')
        a.output.mkdir(parents=True)
        manifest=dict(executable=str(exe),executable_sha256=hashlib.sha256(Path(exe).read_bytes()).hexdigest(),
                      note='Same geometry/strengths. dt_c=.15*dx/c times dt-scale matches; other restrictions may differ.',
                      plans=dict(plans),grid_dx=a.grid_dx)
        (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    failed=[]
    for name,cmd in plans:
        folder=a.output/name
        if not a.analyze_only:
            folder.mkdir()
            (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
            print(f'Running {name}; log: {folder/"run.log"}',flush=True)
            start=time.monotonic()
            with (folder/'run.log').open('w') as log:
                proc=subprocess.run(cmd,cwd=folder,stdout=log,stderr=subprocess.STDOUT)
            (folder/'execution.json').write_text(json.dumps(dict(returncode=proc.returncode,
                elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
            if proc.returncode:
                print(f'{name} failed; preserving partial output',flush=True);failed.append(name);continue
        try:
            result=inspect(folder,a.grid_dx)
            print(json.dumps(dict(case=name,interface=result['interface'],
                                  max_density_error=result['max_density_error'],
                                  max_mach=result['max_mach']),indent=2),flush=True)
        except (OSError,ValueError,RuntimeError,KeyError) as error:
            failed.append(name)
            (folder/'analysis_error.json').write_text(json.dumps(dict(error=str(error)))+'\n')
            print(f'{name}: analysis failed; raw results preserved: {error}',flush=True)
    # analyze-only deliberately refreshes the diagnostic archive; simulations are never overwritten.
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(a.output.rglob('*')):
            if path.is_file() and path.suffix in ('.csv','.json','.log') and 'output' not in path.relative_to(a.output).parts:
                z.write(path,str(Path(a.output.name)/path.relative_to(a.output)))
    print(f'Send back: {archive}',flush=True)
    if failed:raise SystemExit(1)

if __name__=='__main__':main()
