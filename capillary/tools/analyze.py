#!/usr/bin/env python3
"""Independent shape/pressure/load diagnostics; execution success is not physics validation."""
import argparse,csv,json,math,pathlib
import numpy as np
from young_laplace import bridge_reference
from force_audit import force_audit
from scipy.spatial import ConvexHull, QhullError
from scipy.optimize import least_squares

def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def sphere_fit(rows,dx,kind):
 pts=np.array([[float(r[k]) for k in ['x','y','z']] for r in rows if int(r['surface'])==1])
 if kind=='sessile':pts=pts[pts[:,2]>2*dx]
 if len(pts)<12:return {'status':'insufficient independent surface points'}
 A=np.c_[2*pts,np.ones(len(pts))];b=(pts*pts).sum(axis=1)
 x,_,rank,_=np.linalg.lstsq(A,b,rcond=None)
 if rank<4:return {'status':'rank deficient'}
 c=x[:3];R=math.sqrt(max(0,x[3]+c@c));rms=float(np.sqrt(np.mean((np.linalg.norm(pts-c,axis=1)-R)**2)))
 # Surface particle centers lie inside the actual interface. Report both estimates explicitly.
 out={'status':'estimated','points':len(pts),'center':c.tolist(),'radius_particle_centers':R,'radius_half_spacing_offset':R+dx/2,'fit_rms':rms}
 if kind=='sessile':
  out['contact_angle_deg_particle_centers']=math.degrees(math.acos(np.clip(-c[2]/R,-1,1)))
  out['contact_angle_deg_half_spacing_offset']=math.degrees(math.acos(np.clip(-c[2]/(R+dx/2),-1,1)))
 return out


def exterior_sphere_fit(rows, dx, kind):
 """Convex envelope of particle centers, not the library's thick surface band.

 Valid only for a single, approximately convex drop. Offsets are sampling
 sensitivity tests, not a calibrated physical interface reconstruction.
 """
 pts = np.array([[float(r[k]) for k in ('x','y','z')] for r in rows])
 if len(pts) < 12 or not np.isfinite(pts).all():
  return {'status': 'insufficient or invalid points'}
 try:
  hull = pts[ConvexHull(pts).vertices]
 except QhullError:
  return {'status': 'degenerate convex envelope'}
 fits = []
 for cutoff in ([1.5, 2., 2.5] if kind == 'sessile' else [0.]):
  surf = hull[hull[:,2] > cutoff*dx] if kind == 'sessile' else hull
  if len(surf) < 12: continue
  A = np.c_[2*surf, np.ones(len(surf))]
  x, _, rank, _ = np.linalg.lstsq(A, (surf*surf).sum(axis=1), rcond=None)
  if rank < 4: continue
  radius2 = x[3] + x[:3]@x[:3]
  if radius2 <= 0: continue
  fit = least_squares(lambda v: np.linalg.norm(surf-v[:3],axis=1)-v[3],
                      np.r_[x[:3], math.sqrt(radius2)])
  c, radius = fit.x[:3], float(fit.x[3])
  if not fit.success or radius <= 0: continue
  rms = float(np.sqrt(np.mean(fit.fun**2)))
  entry = {'wall_cutoff_dx': cutoff, 'points': len(surf), 'center': c.tolist(),
           'radius_particle_centers': radius, 'fit_rms': rms}
  if kind == 'sessile':
   if abs(c[2]) >= radius: continue
   entry['angle_deg'] = math.degrees(math.acos(-c[2]/radius))
   entry['angle_offset_sensitivity_deg'] = [math.degrees(math.acos(np.clip(-c[2]/(radius+d*dx),-1,1))) for d in (0,.5,1)]
  fits.append(entry)
 if not fits: return {'status': 'no resolved cap fit'}
 primary = min(fits, key=lambda f: abs(f['wall_cutoff_dx']-2))
 out = dict(primary, status='estimated', method='convex envelope geometric sphere fit', cutoff_sensitivity=fits,
            limitation='Assumes one approximately convex drop; offsets are sensitivity, not an error bound. Inspect shape before interpreting.')
 if primary['fit_rms'] > .5*dx: out['status'] = 'poor spherical fit'
 return out

def angle_history(folder, meta):
 index = folder/'liquid_snapshots.csv'
 if not index.exists() or meta['case'] != 'sessile': return None
 records = []
 for row in read(index):
  fit = exterior_sphere_fit(read(folder/row['file']),meta['dx'],'sessile')
  records.append({'step': int(row['step']), 'time': float(row['time']),
    'angle_deg': fit.get('angle_deg'), 'fit_rms': fit.get('fit_rms'), 'status': fit['status']})
 if not records: return None
 with (folder/'angle_history.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
 # Report observed change over the final quarter, without declaring equilibrium.
 tail=[r for r in records if r['time'] >= .75*records[-1]['time'] and r['status']=='estimated']
 out={'samples':len(records),'last':records[-1], 'equilibrium_established':False}
 if len(tail)>1:
  times=np.array([r['time'] for r in tail]); angles=np.array([r['angle_deg'] for r in tail])
  out.update(tail_samples=len(tail),tail_angle_range_deg=float(np.ptp(angles)),
             tail_slope_deg_per_time=float(np.polyfit(times-times[0],angles,1)[0]))
 return out

def analyze(folder):
 folder=pathlib.Path(folder);meta=json.loads((folder/'run.json').read_text());hist=read(folder/'history.csv');last=hist[-1];first=hist[0]
 def vec(r,keys):return np.array([float(r[k]) for k in keys])
 pm=vec(first,['px','py','pz'])+vec(first,['support_impulse_x','support_impulse_y','support_impulse_z'])
 lm=vec(first,['lx','ly','lz'])+vec(first,['support_angular_impulse_x','support_angular_impulse_y','support_angular_impulse_z'])
 p_errors=[];l_errors=[]
 for r in hist:
  p_errors.append(np.linalg.norm(vec(r,['px','py','pz'])+vec(r,['support_impulse_x','support_impulse_y','support_impulse_z'])-pm))
  l_errors.append(np.linalg.norm(vec(r,['lx','ly','lz'])+vec(r,['support_angular_impulse_x','support_angular_impulse_y','support_angular_impulse_z'])-lm))
 out={'metadata':meta,'execution_completed':meta['completed_time'],'physics_status':'NOT_VALIDATED',
 'mass_relative_error':float(last['liquid_mass'])/float(first['liquid_mass'])-1,
 'volume_relative_change':float(last['liquid_volume'])/float(first['liquid_volume'])-1,
 'max_mach':max(float(r['mach']) for r in hist),'max_density_error':max(float(r['max_liquid_density_error']) for r in hist),
 'support_corrected_linear_drift':float(max(p_errors)),'support_corrected_angular_drift':float(max(l_errors)),
 'load_sampling_note':'Upstream wall forces sampled before fluid substep; drift may include temporal splitting and unmatched wetting reactions; the cause is not isolated. No global correction.'}
 snaps=sorted((p for p in folder.glob('liquid_*.csv') if p.stem.split('_')[-1].isdigit()),key=lambda p:int(p.stem.split('_')[-1]));rows=read(snaps[-1])
 if meta['case'] in ('laplace','sessile'):
  out['legacy_surface_band_fit']=sphere_fit(rows,meta['dx'],meta['case'])
  shape=exterior_sphere_fit(rows,meta['dx'],meta['case']);out['shape']=shape
  if shape['status']=='estimated':
   if meta['case']=='laplace':
    dp=float(last['pressure_inside'])-float(last['pressure_outside']);ref=2*meta['sigma']/(shape['radius_particle_centers']+meta['dx']/2)
    out['laplace']={'measured':dp,'reference':ref,'relative_error':(dp-ref)/ref}
   else:out['angle_error_deg']=shape['angle_deg']-meta['theta']
 elif meta['case']=='bridge':
  try:
   ref=bridge_reference(float(last['liquid_volume']),meta['gap'],meta['theta'],meta['radius'],meta['sigma']);out['young_laplace']=ref
   loads=read(folder/'loads.csv');top=[float(r['fz']) for r in loads if r['body']=='1'];tail=top[-max(1,len(top)//5):]
   if tail:out['bridge_force']={'mean_top_fz':float(np.mean(tail)),'top_fz_std':float(np.std(tail)),'reference_top_fz':-ref['attractive_force'],'relative_error':float((np.mean(tail)+ref['attractive_force'])/ref['attractive_force'])}
  except (RuntimeError,ValueError) as e:out['young_laplace']={'status':'reference_not_resolved','reason':str(e)}
 if 'liquid_max_speed' in first:
  out['max_liquid_mach']=max(float(r['liquid_max_speed'])/meta['sound_speed'] for r in hist)
  out['max_gas_mach']=max(float(r['gas_max_speed'])/meta['sound_speed'] for r in hist)
  out['max_gas_density_error']=max(float(r['gas_density_error']) for r in hist)
 out['angle_history']=angle_history(folder,meta)
 out['force_audit']=force_audit(folder)
 (folder/'analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('directory',type=pathlib.Path);a=p.parse_args();print(json.dumps(analyze(a.directory),indent=2))
