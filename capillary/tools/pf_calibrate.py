#!/usr/bin/env python3
"""Stage-gated PF-SPH calibration: multi-radius Laplace -> sessile angle -> two-dx bridge.

Runs the native executable on the user's computer. No external data is fetched.
Existing results are never overwritten. Standalone --dry-run prints exact commands.
"""
import argparse
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
import zipfile
from pf_calibration_core import (analyze_droplet,fit_surface_tension,
    evaluate_angle_run,select_angle,bridge_summary,compare_bridges)


def positive(v):
    x=float(v)
    if not (math.isfinite(x) and x>0):raise argparse.ArgumentTypeError('positive finite number required')
    return x


def nonnegative(v):
    x=float(v)
    if not (math.isfinite(x) and x>=0):raise argparse.ArgumentTypeError('finite nonnegative number required')
    return x


def exe_path(value,dry):
    path=shutil.which(value)
    if not path and Path(value).is_file():path=str(Path(value).resolve())
    if not path and not dry:raise ValueError(f'Cannot locate executable: {value}')
    return path or value


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def read(path):
    return json.loads(Path(path).read_text())


def launch(cmd,folder,dry=False):
    print(shlex.join(list(map(str,cmd))),flush=True)
    if dry:return 0
    folder.mkdir(parents=True,exist_ok=False)
    write(folder/'command.json',cmd)
    start=time.monotonic()
    with (folder/'run.log').open('w') as f:
        rc=subprocess.run(cmd,cwd=folder,stdout=f,stderr=subprocess.STDOUT).returncode
    write(folder/'execution.json',dict(exit_code=rc,wall_seconds=time.monotonic()-start))
    if rc:
        print(f'Native run exited with {rc}: {folder}. Raw results preserved.',file=sys.stderr)
    return rc


def archive_output(root):
    """Compact report archive. Preserve full checkpoint/snapshots locally for restart."""
    root=Path(root).resolve();dest=root.with_suffix('.zip')
    if dest.exists():
        print(f'Archive already exists and is preserved: {dest}. New analysis remains in the directory.',flush=True)
        return
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(root.rglob('*')):
            if p.is_file() and p.suffix in ('.json','.csv','.log','.md') and not p.name.startswith('pf_snapshot_'):
                z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
    print('Send diagnostic archive:',dest,flush=True)


def common(p):
    p.add_argument('--executable',default='build/pf_capillary_validation')
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--dx',type=positive,required=True)
    p.add_argument('--end-time',type=positive,default=10.)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--dt-scale',type=positive,default=.75)
    p.add_argument('--record-every',type=int,default=50)
    p.add_argument('--snapshot-interval',type=nonnegative,default=1.)
    p.add_argument('--max-steps',type=int,default=1000000)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--analyze-only',action='store_true')


def base_cmd(a,case,radius,sff,sfs):
    return [exe_path(a.executable,a.dry_run),'--case',case,'--dx',str(a.dx),'--radius',str(radius),
            '--sff',str(sff),'--sfs',str(sfs),'--sound-speed','10','--viscosity','.05',
            '--dt-scale',str(a.dt_scale),'--end-time',str(a.end_time),'--steps',str(a.max_steps),
            '--threads',str(a.threads),'--record-every',str(a.record_every),
            '--dump-every','0','--snapshot-interval',str(a.snapshot_interval),
            '--checkpoint-interval',str(min(1.,a.end_time))]


def check_args(a):
    if a.threads<=0 or a.record_every<=0 or a.max_steps<=0 or a.dt_scale>1:
        raise ValueError('Threads, record rate, and steps must be positive and dt-scale <=1')
    if not a.analyze_only and not a.dry_run and (a.output.exists() or a.output.with_suffix('.zip').exists()):
        raise ValueError('Output already exists: choose a fresh --output or use --analyze-only')
    if a.analyze_only and not a.output.is_dir():raise ValueError('Analyze-only requires existing output directory')


def run_droplets(a):
    check_args(a)
    if len(a.radii)<3 or len(set(a.radii))!=len(a.radii):
        raise ValueError('Provide at least three distinct initial droplet radii')
    if min(a.radii)<=12*a.dx:  # ensures R-2H>4dx initially
        raise ValueError('Use R > 12 dx to resolve a bulk virial-pressure core')
    plans=[(f'droplet_r{r:g}'.replace('.','p'),base_cmd(a,'droplet',r,a.sff,0.)) for r in a.radii]
    if a.dry_run:
        for name,cmd in plans: print(name,shlex.join(cmd))
        return
    if not a.analyze_only:
        a.output.mkdir(parents=True)
        write(a.output/'manifest.json',dict(stage='droplet',dx=a.dx,sff=a.sff,radii=a.radii,
             note='Measure total virial pressure vs 2/volume-equivalent radius; do not use paper mapping as data.'))
    reports=[];failed=[]
    for name,cmd in plans:
        folder=a.output/name
        if not a.analyze_only and launch(cmd,folder):failed.append(name);continue
        try:
            summary=analyze_droplet(folder)
            write(folder/'droplet_analysis.json',summary)
            reports.append(summary)
        except Exception as exc:
            failed.append(f'{name}: {exc}')
    if len(reports)>=3:
        fit=fit_surface_tension(reports)
        if a.target_gamma is not None:
            error=abs(fit['gamma']-a.target_gamma)/a.target_gamma
            fit['target_gamma']=a.target_gamma
            fit['target_gamma_relative_error']=error
            fit['suggested_next_sff']=a.sff*a.target_gamma/fit['gamma'] if fit['gamma']>0 else None
            if error>a.gamma_tolerance and fit['status']=='provisional_pass':
                fit['status']='target_mismatch'
                fit['note']='Repeat droplet calibration at suggested_next_sff; do not reuse previous gamma as measured.'
        write(a.output/'gamma_report.json',fit)
        print(json.dumps({k:v for k,v in fit.items() if k!='droplets'},indent=2))
    else:
        write(a.output/'gamma_report.json',dict(status='insufficient_droplets',failure=failed,droplets=reports))
    archive_output(a.output)
    if failed:
        raise ValueError(f'Droplet cases incomplete: {failed}. Inspect logs or continue runs.')


def run_angles(a):
    check_args(a)
    gamma=read(a.gamma_report)
    if gamma.get('status')!='provisional_pass':
        raise ValueError('Gamma report not yet converged; inspect it before spending compute on angles')
    if not math.isclose(gamma['dx'],a.dx,abs_tol=1e-12):
        raise ValueError('Gamma report and sessile dx do not match')
    sff=gamma['sff']
    if a.sfs is not None and a.sfs_ratios is not None:
        raise ValueError('Specify either --sfs or --sfs-ratios, not both')
    sfs=a.sfs if a.sfs is not None else [r*sff for r in
          (a.sfs_ratios if a.sfs_ratios is not None else [1.,1.5,2.,3.,4.])]
    if not sfs or any(v<=0 for v in sfs) or len(set(sfs))!=len(sfs):
        raise ValueError('Supply distinct positive fluid-solid interaction strengths')
    plans=[(f'sessile_sfs{s:.8g}'.replace('.','p'),base_cmd(a,'sessile',a.drop_radius,sff,s)) for s in sfs]
    if a.dry_run:
        for name,cmd in plans:print(name,shlex.join(cmd))
        return
    if not a.analyze_only:
        a.output.mkdir(parents=True)
        write(a.output/'manifest.json',dict(stage='angle',dx=a.dx,sff=sff,trial_sfs=sfs,
            target_angle_deg=a.target_angle,gamma_report=str(a.gamma_report)))
    reports=[];failed=[]
    for name,cmd in plans:
        folder=a.output/name
        if not a.analyze_only and launch(cmd,folder): failed.append(name);continue
        try:
            report=evaluate_angle_run(folder,a.grid_dx)
            reports.append(report)
            write(folder/'angle_analysis.json',report)
        except Exception as exc:
            failed.append(f'{name}: {exc}')
    write(a.output/'angle_report.json',dict(target_angle_deg=a.target_angle,dx=a.dx,
        reports=reports,failures=failed))
    try:
        selected=select_angle(reports,gamma,a.target_angle,a.angle_tolerance)
        selected['gamma_report_source']=str(a.gamma_report)
        write(a.output/'calibration.json',selected)
        print('Provisional coefficients:',json.dumps(selected,indent=2))
    except ValueError as exc:
        write(a.output/'needs_refinement.json',dict(reason=str(exc),target_angle_deg=a.target_angle,
              note='Try extra sfs values or extend relaxation; never silently pick an unstable angle.'))
        print('Angle selection unresolved:',exc,flush=True)
    archive_output(a.output)
    if failed:raise ValueError(f'Sessile runs unresolved: {failed}')


def run_bridge(a):
    check_args(a)
    cal=read(a.calibration)
    if cal.get('status')!='provisional_pass' or not math.isclose(cal['dx'],a.dx,abs_tol=1e-12):
        raise ValueError('Provide converged provisional calibration for this exact dx')
    cmd=base_cmd(a,'bridge',a.liquid_radius,cal['sff'],cal['sfs'])
    cmd+=['--sphere-radius',str(a.sphere_radius),'--bridge-gap',str(a.gap),
          '--fixed-spheres','1','--contact-stiffness','0']
    if a.dry_run:print(shlex.join(cmd));return
    a.output.mkdir(parents=True)
    write(a.output/'calibration_used.json',cal)
    if launch(cmd,a.output/'bridge'):
        raise ValueError('Native bridge run failed; inspect the preserved results')
    summary=bridge_summary(a.output/'bridge')
    write(a.output/'bridge_analysis.json',summary)
    archive_output(a.output)
    print(json.dumps(summary,indent=2))


def run_compare(a):
    coarse=bridge_summary(a.coarse)
    fine=bridge_summary(a.fine)
    report=compare_bridges(coarse,fine,read(a.coarse_cal),read(a.fine_cal),a.gamma_tol,a.angle_tol,a.mass_tol)
    report.update(coarse=coarse,fine=fine)
    write(a.output,report)
    print(json.dumps(report,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='stage',required=True)
    droplet=sub.add_parser('droplets',help='Measure gamma from >=3 equilibrated free droplets')
    common(droplet)
    droplet.add_argument('--radii',nargs='+',type=positive,default=[.65,.8,.95])
    droplet.add_argument('--sff',type=positive,default=.001)
    droplet.add_argument('--target-gamma',type=positive,default=None,
         help='For refinement: previously measured target gamma at coarser dx')
    droplet.add_argument('--gamma-tolerance',type=positive,default=.05)
    wet=sub.add_parser('angles',help='Measure theta(sfs/sff) with stationary sessile droplets')
    common(wet)
    wet.add_argument('--gamma-report',required=True,type=Path)
    wet.add_argument('--target-angle',type=positive,default=60.)
    wet.add_argument('--angle-tolerance',type=positive,default=5.)
    wet.add_argument('--sfs',nargs='+',type=positive,default=None)
    wet.add_argument('--sfs-ratios',nargs='+',type=positive,default=None,
         help='Ratios sfs/sff. Default: 1 1.5 2 3 4')
    wet.add_argument('--drop-radius',type=positive,default=.6)
    wet.add_argument('--grid-dx',type=positive,default=.3)
    bridge=sub.add_parser('bridge',help='Run fixed bridge only with calibrated coefficients')
    common(bridge)
    bridge.add_argument('--calibration',type=Path,required=True)
    bridge.add_argument('--sphere-radius',type=positive,default=.3)
    bridge.add_argument('--gap',type=positive,default=.2)
    bridge.add_argument('--liquid-radius',type=positive,default=.35)
    compare=sub.add_parser('compare',help='Compare two independently calibrated bridge resolutions')
    compare.add_argument('--coarse',required=True,type=Path)
    compare.add_argument('--fine',required=True,type=Path)
    compare.add_argument('--coarse-cal',required=True,type=Path)
    compare.add_argument('--fine-cal',required=True,type=Path)
    compare.add_argument('--output',required=True,type=Path)
    compare.add_argument('--gamma-tol',type=positive,default=.1)
    compare.add_argument('--angle-tol',type=positive,default=5.)
    compare.add_argument('--mass-tol',type=positive,default=.04)
    a=p.parse_args()
    try:
        if a.stage=='droplets':run_droplets(a)
        elif a.stage=='angles':run_angles(a)
        elif a.stage=='bridge':run_bridge(a)
        else:run_compare(a)
    except (OSError,ValueError,RuntimeError) as exc:
        p.exit(1,f'Calibration {a.stage} failed: {exc}\n')


if __name__=='__main__':main()
