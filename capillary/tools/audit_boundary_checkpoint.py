"""Read-only PFDEM checkpoint-prefix audit; does not restore/evaluate DEM springs.
Matches pf_liquid::accumulate_pair, using the stored liquid/rigid boundary state.
Requires NumPy. Reports the saved checkpoint time, not the last failed log time.
"""
import argparse
import numpy as np,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('checkpoint',type=Path)
parser.add_argument('--output',type=Path)
args=parser.parse_args()
lines=args.checkpoint.read_text().splitlines()
if not lines or lines[0] not in ('PFDEM 1','PFDEM 2','PFDEM 3','PFDEM 4','PFDEM 5'):raise ValueError('Unsupported checkpoint format')
p=np.array(list(map(float,lines[2].split())));dx,rho0,c,mu,sff,sfs=p;H=4*dx
header=lines[4].split();nf,n,nb,ns,time,step=header[:6];nf,n,nb,ns,step=map(int,(nf,n,nb,ns,step));start=11
b=np.array([list(map(float,x.split())) for x in lines[start:start+nb]]);start+=nb
s=np.array([list(map(float,x.split())) for x in lines[start:start+n]]);start+=n
sites=np.array([list(map(float,x.split())) for x in lines[start:start+ns]])
X=s[:,:3];V=s[:,3:6];rho=s[:,6];mass=s[:,7];solid=s[:,8].astype(bool)
owner=np.full(n,-1,int);owner[sites[:,0].astype(int)]=sites[:,1].astype(int)
exposed=np.ones(n,bool)
weights=np.ones(n);partition=lines[0]=='PFDEM 3' or (lines[0] in ('PFDEM 4','PFDEM 5') and header[6]=='0')
for site in sites:
 idx,own=map(int,site[:2])
 for j in range(nb):
  if j==own:continue
  if np.linalg.norm(b[j,:3]-b[own,:3])>=b[j,13]+b[own,13]:continue
  if np.linalg.norm(X[idx]-b[j,:3])<b[j,13]:
   if partition:weights[idx]+=1
   else:exposed[idx]=False;break
 if partition:weights[idx]=1./weights[idx]
P=c*c*rho0/7*((rho/rho0)**7-1)
gaps=np.linalg.norm(X[:nf,None,:]-b[None,:,:3],axis=2)-b[None,:,13]
allgaps=gaps.min(axis=1);report={'step':step,'time_s':float(time),'fluid_count':nf,'rho_ratio_percentiles':np.percentile(rho[:nf]/rho0,[0,25,50,75,100]).tolist(),'inside_fluid_count':int((allgaps<0).sum()),'deeper_than_dx_count':int((allgaps < -dx).sum()),'min_surface_gap_nm':float(allgaps.min()*1e9),'parameters':dict(zip(['dx','rho0','sound','mu','sff','sfs'],p)),'probes':[],'overlap_rule':('volume_partition' if partition else 'legacy_mask')}
for i in [int(np.argmin(rho[:nf])),int(np.argmax(rho[:nf]))]:
 d=X[i]-X;r=np.linalg.norm(d,axis=1)
 distinct=np.arange(n)!=i
 if np.any((r==0)&distinct&exposed):raise ValueError(f'Coincident active SPH samples at probe {i}; solver rejects this state too')
 active=distinct&(r<H)&exposed
 idx=np.flatnonzero(active);dist=r[idx];e=d[idx]/dist[:,None];q=dist/H;z=1-q
 dw=-210/(np.pi*H**4)*q*z**3;radial=((V[i]-V[idx])*e).sum(axis=1)
 vol=mass[i]*mass[idx]/(rho[i]*rho[idx]);a=np.where(solid[idx],-11.,-23.5)
 profile=z*z*(1+2*q+(3+a)*q*q)
 components={'pressure':(-vol*(P[i]+P[idx])*dw)[:,None]*e,'pair_force':(H*np.where(solid[idx],sfs,sff)*profile)[:,None]*e,'viscous':(10*mu*vol*dw/dist*radial)[:,None]*e}
 components={name:force*weights[idx,None] for name,force in components.items()}
 body=int(np.argmin(gaps[i]));normal=(X[i]-b[body,:3]);normal_length=np.linalg.norm(normal)
 normal_defined=bool(normal_length>0)
 if normal_defined:normal/=normal_length
 # A sample exactly at the sphere center has no unique outward normal.
 def project(v):return float(v@normal) if normal_defined else None
 inertia=.4*b[body,14]*b[body,13]**2;omega=b[body,6:9]/inertia
 u=V[i]-b[body,3:6]-np.cross(omega,X[i]-b[body,:3])
 out={'sample':i,'rho_ratio':rho[i]/rho0,'nearest_body':body,'gap_nm':gaps[i,body]*1e9,'relative_outward_velocity_m_s':project(u),'nearest_boundary_sample_distance_dx':(float(dist[solid[idx]].min()/dx) if solid[idx].any() else None),'components':{}}
 for group,mask in [('liquid',~solid[idx]),('solid',solid[idx]),('nearest_body_sites',owner[idx]==body)]:
  sub={}
  for name,force in components.items():sub[name]={'force_N':force[mask].sum(axis=0).tolist(),'outward_force_N':project(force[mask].sum(axis=0))}
  rate=(weights[idx]*rho[i]*mass[idx]/rho[idx]*radial*dw)[mask].sum()
  sub['continuity_divergence_s_inv']=float(-rate/rho[i]);sub['count']=int(mask.sum());out['components'][group]=sub
 out['total_force_N']=sum(v.sum(axis=0) for v in components.values()).tolist()
 # Same active edges and per-pair formulas as solver; no EOS clipping.
 report['probes'].append(out)
if args.output:
 if args.output.exists():raise ValueError('Audit output exists; use a new path')
 args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(f"checkpoint step={step} time_ns={float(time)*1e9:.8g} liquid={nf} inside={report['inside_fluid_count']} deeper_than_dx={report['deeper_than_dx_count']} min_gap_nm={report['min_surface_gap_nm']:.8g}")
for probe in report['probes']:
 print(f"sample={probe['sample']} rho/rho0={probe['rho_ratio']:.8g} gap_nm={probe['gap_nm']:.8g} body={probe['nearest_body']} outward_relative_speed={probe['relative_outward_velocity_m_s']}")
 for group in ('liquid','solid'):
  values=probe['components'][group]
  print(group, 'outward force [pN]:', {k:(round(v['outward_force_N']*1e12,8) if v['outward_force_N'] is not None else None) for k,v in values.items() if isinstance(v,dict)}, 'divergence [1/s]:',values['continuity_divergence_s_inv'])
