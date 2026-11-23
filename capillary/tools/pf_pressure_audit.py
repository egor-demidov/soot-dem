#!/usr/bin/env python3
"""Offline virial-pressure audit of saved PF-SPH droplets. Runs no fluid dynamics.

Reports the paper-style per-particle half virial, separating EOS pressure and PF
forces. Core-radius curves and quadrature-volume normalization are diagnostics;
none is substituted for the production calibration or used to force agreement.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def particle_virials(xyz,rho,meta,chunk_size=65536):
    xyz=np.asarray(xyz,dtype=float);rho=np.asarray(rho,dtype=float)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or rho.shape!=(len(xyz),) or not len(xyz):
        raise ValueError('Invalid particle arrays')
    if not np.isfinite(xyz).all() or not np.isfinite(rho).all() or np.any(rho<=0):
        raise ValueError('Nonfinite state or nonpositive density')
    if chunk_size<1:raise ValueError('Positive chunk size required')
    dx=float(meta['dx']);rho0=float(meta.get('rho0',1.));c=float(meta['sound_speed']);sff=float(meta['sff'])
    if not all(math.isfinite(v) and v>0 for v in (dx,rho0,c)) or not(math.isfinite(sff) and sff>=0):
        raise ValueError('Invalid PF material parameters')
    H=4*dx;mass=rho0*dx**3
    pressure=c*c*rho0/7*((rho/rho0)**7-1)
    pairs=cKDTree(xyz).query_pairs(H,output_type='ndarray')
    eos=np.zeros(len(xyz));pf=np.zeros(len(xyz))
    for begin in range(0,len(pairs),chunk_size):
        ids=pairs[begin:begin+chunk_size];i,j=ids[:,0],ids[:,1]
        d=xyz[i]-xyz[j];r=np.sqrt(np.einsum('ij,ij->i',d,d))
        if np.any(r==0):raise ValueError('Coincident distinct fluid samples')
        q=r/H;z=np.maximum(0,1-q)
        dw=-210/(math.pi*H**4)*q*z**3
        fp=-mass*mass/(rho[i]*rho[j])*(pressure[i]+pressure[j])*dw
        fc=H*sff*z*z*(1+2*q-20.5*q*q)
        for dest,force in ((eos,fp),(pf,fc)):
            virial=.5*r*force
            dest+=np.bincount(i,weights=virial,minlength=len(xyz))
            dest+=np.bincount(j,weights=virial,minlength=len(xyz))
    return eos,pf,pressure,mass,len(pairs)


def audit_frame(folder):
    folder=Path(folder);meta=json.loads((folder/'pf_run.json').read_text())
    if meta['case']!='droplet':raise ValueError('Requires isolated liquid droplet, with no solid samples')
    data=np.genfromtxt(folder/'pf_final_liquid.csv',delimiter=',',names=True)
    xyz=np.column_stack([data[k] for k in ('x','y','z')]);rho=np.atleast_1d(data['rho'])
    eos,pf,pressure,mass,pairs=particle_virials(xyz,rho,meta)
    vol=mass/rho;R=float(np.cbrt(3*vol.sum()/(4*math.pi)));H=4*meta['dx']
    center=xyz.mean(axis=0);radii=np.linalg.norm(xyz-center,axis=1)
    rv=R-2*H
    if rv<=0:raise ValueError('Unresolved R-2H core')
    scan=np.unique(np.r_[np.linspace(max(.25*H,.1*R),R-H,32),rv])
    records=[]
    for radius in scan:
        inside=radii<radius;n=int(inside.sum())
        if not n:continue
        V=4*math.pi*radius**3/3;Vq=float(vol[inside].sum())
        pe=float(eos[inside].sum()/(3*V));pc=float(pf[inside].sum()/(3*V))
        records.append(dict(core_radius=float(radius),core_radius_over_H=float(radius/H),
            interface_margin_over_H=float((R-radius)/H),particles=n,
            geometric_volume=V,quadrature_volume=Vq,quadrature_volume_relative_error=Vq/V-1,
            pressure_eos_pair=pe,pressure_pf_pair=pc,pressure_total=pe+pc,
            pressure_total_quadrature_normalized=(pe+pc)*V/Vq,
            mean_eos_pressure=float(pressure[inside].mean()),
            cancellation_factor=(abs(pe)+abs(pc))/max(abs(pe+pc),1e-30)))
    fixed=min(records,key=lambda r:abs(r['core_radius']-rv))
    with (folder/'pf_history.csv').open() as f:history=list(csv.DictReader(f))
    last=history[-1]
    stored=float(last['virial_pressure'])
    # The CSV frame and history endpoint must refer to the same completed state.
    matching=bool(meta.get('completed_time') and math.isclose(float(last['time']),float(meta['time']),rel_tol=1e-10,abs_tol=1e-12))
    band=[r['pressure_total'] for r in records if R>4*H and 2-1e-10<=r['interface_margin_over_H']<=3+1e-10]
    result=dict(case=folder.name,dx=meta['dx'],sff=meta['sff'],volume_radius=R,support=H,
        time=meta['time'],radius_over_support=R/H,pair_count=pairs,
        frame_matches_history_endpoint=matching,fixed_core=fixed,
        stored_final_pressure=stored,independent_pressure_difference=fixed['pressure_total']-stored,
        radial_band_mean=float(np.mean(band)) if band else None,
        radial_band_range=float(np.ptp(band)) if band else None,
        warning='One saved frame is not a time-averaged calibration. Alternative volume normalization is diagnostic only.',
        radial_profile=records)
    return result


def audit_directory(root,output,plot=False):
    root=Path(root);output=Path(output)
    folders=sorted(p.parent for p in root.glob('*/pf_run.json'))
    if not folders:raise ValueError('No droplet case directories found')
    reports=[audit_frame(folder) for folder in folders]
    output.mkdir(parents=True,exist_ok=True)
    for report in reports:
        (output/(report['case']+'.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    result=dict(status='diagnostic_only',physics_validated=False,cases=reports,
        note='Independent final-frame pressure reconstruction and spatial-sampling audit; run no trajectories.')
    (output/'pressure_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        for report in reports:
            rows=report['radial_profile'];label=f"R/H={report['radius_over_support']:.2f}"
            axes[0].plot([r['interface_margin_over_H'] for r in rows],[r['pressure_total'] for r in rows],label=label)
            axes[1].plot([r['interface_margin_over_H'] for r in rows],[100*r['quadrature_volume_relative_error'] for r in rows],label=label)
        for ax in axes:
            ax.axvspan(2,3,color='grey',alpha=.12);ax.set_xlabel('Distance from volume radius / H');ax.legend()
        axes[0].set_ylabel('Conservative virial pressure (simulation units)')
        axes[1].set_ylabel('Sample-volume / geometric-volume error (%)')
        fig.savefig(output/'pressure_core_sensitivity.png',dpi=160);plt.close(fig)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--plot',action='store_true')
    a=p.parse_args();result=audit_directory(a.directory,a.output,a.plot)
    print(json.dumps([dict(case=r['case'],stored=r['stored_final_pressure'],reconstructed=r['fixed_core']['pressure_total'],
        difference=r['independent_pressure_difference'],radial_band_mean=r['radial_band_mean'],radial_band_range=r['radial_band_range']) for r in result['cases']],indent=2))
