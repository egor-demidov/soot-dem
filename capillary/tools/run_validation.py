#!/usr/bin/env python3
"""Run in fresh case directories, analyze, and archive evidence. Never launches a build."""
import argparse,hashlib,json,pathlib,shutil,subprocess,sys,time,zipfile
from analyze import analyze
ROOT=pathlib.Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('--executable',type=pathlib.Path,default=ROOT/'build/capillary_validation')
p.add_argument('--output',type=pathlib.Path,default=ROOT/'results-validation')
p.add_argument('--smoke',action='store_true',help='Only two acoustic steps per case; no physical validation')
p.add_argument('--refinement',action='store_true',help='Three resolutions, three wetting angles and sound-speed sensitivity')
p.add_argument('--end-time',type=float,default=20.)
p.add_argument('--threads',type=int,default=4)
p.add_argument('--timeout',type=int,default=86400,help='Wall seconds per case')
p.add_argument('--cases',nargs='+',choices=['laplace','sessile','bridge','moving'],default=['laplace','sessile','bridge','moving'])
a=p.parse_args();exe=a.executable.resolve();out=a.output.resolve()
if not exe.is_file():p.error('Build capillary_validation first')
if out.exists():p.error('Output exists: select a fresh --output directory')
if a.threads<1 or a.end_time<=0:p.error('Positive threads and end time required')
out.mkdir(parents=True)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
manifest={'executable_sha256':digest(exe),'sources':{str(f.relative_to(ROOT)):digest(f) for folder in ['src','tools'] for f in (ROOT/folder).rglob('*') if f.is_file() and f.suffix in ['.cpp','.py']},'smoke_only':a.smoke,'cases':[],'pin':'319520e02f631847f85160c0e8b55648c395eafa'}
failed=False
try:
 for kind in a.cases:
  sizes=[.3] if a.smoke else ([.12,.08,.06] if a.refinement else [.12])
  angles=[60,90,120] if kind=='sessile' and a.refinement and not a.smoke else [60]
  for dx in sizes:
   for theta in angles:
    settings=[(10.,1.),(20.,1.),(10.,.5)] if a.refinement and not a.smoke and dx==sizes[-1] else [(10.,1.)]
    for sound,dt_scale in settings:
     name=f'{kind}_dx{dx:g}_theta{theta}_c{sound:g}_dt{dt_scale:g}';dest=out/name;dest.mkdir()
     cmd=[str(exe),'--case',kind,'--dx',str(dx),'--theta',str(theta),'--sound-speed',str(sound),'--dt-scale',str(dt_scale),'--end-time',str(.05 if a.smoke else a.end_time),'--steps',str(2 if a.smoke else 10000000),'--threads',str(a.threads),'--record-every',str(1 if a.smoke else 20),'--dump-every',str(0 if a.smoke else 1000)]
     entry={'name':name,'command':cmd};manifest['cases'].append(entry)
     print('Running',name,flush=True);start=time.monotonic()
     try:
      with (dest/'run.log').open('w') as log:proc=subprocess.run(cmd,cwd=dest,stdout=log,stderr=subprocess.STDOUT,timeout=a.timeout)
      entry['returncode']=proc.returncode
      if proc.returncode==0:entry['analysis']=analyze(dest)
      else:failed=True
     except Exception as e:entry['error']=str(e);failed=True
     entry['wall_seconds']=time.monotonic()-start
     (out/'summary.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
finally:
 (out/'summary.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
 archive=out.with_suffix('.zip')
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
  for f in out.rglob('*'):
   if f.is_file():z.write(f,str(f.relative_to(out.parent)))
 print('Results:',archive)
print('Execution failures:',failed,'; physical validation remains subject to equilibrium and refinement checks.')
sys.exit(1 if failed else 0)
