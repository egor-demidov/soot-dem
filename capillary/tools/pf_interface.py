#!/usr/bin/env python3
"""Offline liquid-only density isosurface and zero-gravity, volume-constrained cap fit.

Whebell et al. (2025), eqs.3.6-3.7. At g=0 the Young-Laplace profile is
an exact spherical cap. This does not implement the finite-gravity ODE fit.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares


def read_samples(path, mass):
    a = np.genfromtxt(path, delimiter=',', names=True)
    xyz = np.column_stack([a[k] for k in ('x','y','z')])
    rho = np.atleast_1d(a['rho'])
    if len(xyz)<12 or not np.isfinite(xyz).all() or not np.isfinite(rho).all() or np.any(rho<=0):
        raise ValueError('Invalid or insufficient liquid samples')
    return xyz, rho, np.full(len(xyz), mass)


def read_vtp_samples(path, mass):
    import xml.etree.ElementTree as ET
    piece=ET.parse(path).find('./PolyData/Piece')
    if piece is None:
        raise ValueError('Expected liquid VTP PolyData')
    positions=piece.find('./Points/DataArray')
    density=piece.find('./PointData/DataArray[@Name="Density"]')
    for array in (positions,density):
        if array is None or array.get('format')!='ascii':
            raise ValueError('Expected ASCII liquid VTP positions and Density')
    count=int(piece.attrib['NumberOfPoints'])
    xyz=np.fromstring(positions.text or '',sep=' ').reshape(-1,3)
    rho=np.fromstring(density.text or '',sep=' ')
    if (count<12 or len(xyz)!=count or len(rho)!=count or
            not np.isfinite(xyz).all() or not np.isfinite(rho).all() or np.any(rho<=0)):
        raise ValueError('Invalid or insufficient liquid samples in VTP')
    mass_array=piece.find('./PointData/DataArray[@Name="Mass"]')
    if mass_array is None:
        masses=np.full(count,mass)
    else:
        if mass_array.get('format')!='ascii':
            raise ValueError('Expected ASCII liquid VTP masses')
        masses=np.fromstring(mass_array.text or '',sep=' ')
        if len(masses)!=count or not np.isfinite(masses).all() or np.any(masses<=0):
            raise ValueError('Invalid liquid masses in VTP')
    return xyz,rho,masses


def density_grid(xyz, masses, support, spacing, max_nodes=30_000_000):
    """Deposit exact compact Wendland kernels onto a grid, without a global neighbor matrix."""
    lo = np.floor((xyz.min(axis=0)-support)/spacing)*spacing
    hi = np.ceil((xyz.max(axis=0)+support)/spacing)*spacing
    shape = np.rint((hi-lo)/spacing).astype(int)+1
    if math.prod(int(v) for v in shape)>max_nodes:
        raise ValueError('Interface grid exceeds memory limit; use larger --grid-dx')
    axes = [lo[k]+spacing*np.arange(shape[k]) for k in range(3)]
    density = np.zeros(tuple(shape), dtype=np.float32)
    factor = 21/(2*math.pi*support**3)
    for point,mass in zip(xyz,masses):
        lower = np.maximum(0,np.ceil((point-support-lo)/spacing).astype(int))
        upper = np.minimum(shape,np.floor((point+support-lo)/spacing).astype(int)+1)
        slices = tuple(slice(int(l),int(u)) for l,u in zip(lower,upper))
        x,y,z = [axes[k][slices[k]]-point[k] for k in range(3)]
        q = np.sqrt(x[:,None,None]**2+y[None,:,None]**2+z[None,None,:]**2)/support
        t = np.maximum(0,1-q)
        density[slices] += (mass*factor*t**4*(1+4*q)).astype(np.float32)
    return density,lo


def reconstruct(xyz,masses,support,spacing,rho0=1.):
    from skimage.measure import marching_cubes
    grid,lo = density_grid(xyz,masses,support,spacing)
    level = rho0/2
    if not grid.min()<level<grid.max(): raise ValueError('No resolved rho0/2 isosurface')
    vertices,faces,_,_ = marching_cubes(grid,level=level,spacing=(spacing,)*3,allow_degenerate=False)
    vertices += lo
    # Center coordinates before volume calculation to reduce translation cancellation.
    v = vertices-vertices.mean(axis=0)
    volume = abs(np.einsum('ij,ij->i',v[faces[:,0]],np.cross(v[faces[:,1]],v[faces[:,2]])).sum())/6
    return vertices,faces,float(volume)


def fit_cap(vertices,volume,dx,center_xy,wall_z=0.,fit_shift=True):
    """Known fluid volume constrains curvature/angle; wall/interface shift is fitted.

    Returns cutoff sensitivity, without selecting a desired target angle or adjusting
    the physical wall to force agreement. Assumes one convex cap at zero gravity.
    """
    if not volume>0: raise ValueError('Positive cap volume required')
    radial = np.linalg.norm(vertices[:,:2]-center_xy,axis=1)
    fits=[]
    for cutoff in (1.,1.5,2.):
        keep=vertices[:,2]>wall_z+cutoff*dx
        r,z=radial[keep],vertices[keep,2]
        if len(r)<30: continue
        # Paper's inverse-radius cylindrical weighting, regularized on the axis.
        weight=np.sqrt(1/(r+.1*dx))
        weight/=np.sqrt(np.mean(weight**2))
        def geometry(x):
            theta=x[0]
            shift=x[1] if fit_shift else 0.
            c=np.cos(theta)
            R=(3*volume/(np.pi*(2-3*c+c**3)))**(1/3)
            return R,wall_z+shift-R*c
        def residual(x):
            R,cz=geometry(x)
            return weight*(np.hypot(r,z-cz)-R)
        # Multiple starts avoid declaring a local minimum as a unique cap fit.
        best=None
        for angle in (40.,90.,140.):
            initial=[np.deg2rad(angle)]
            lower,upper=[np.deg2rad(5)],[np.deg2rad(175)]
            if fit_shift:
                initial.append(0.)
                lower.append(-4*dx)
                upper.append(4*dx)
            opt=least_squares(residual,initial,bounds=(lower,upper),
                    xtol=1e-10,ftol=1e-10,gtol=1e-10)
            if best is None or np.linalg.norm(opt.fun)<np.linalg.norm(best.fun): best=opt
        R,cz=geometry(best.x)
        rms=float(np.sqrt(np.mean((np.hypot(r,z-cz)-R)**2)))
        bound=bool(np.any(best.active_mask))
        fits.append(dict(cutoff_dx=cutoff,vertices=int(len(r)),angle_deg=float(np.rad2deg(best.x[0])),
                    interface_shift=float(best.x[1]) if fit_shift else 0.,radius=float(R),center_z=float(cz),fit_rms=rms,
                    optimizer_success=bool(best.success),at_parameter_bound=bound))
    if not fits: return dict(status='insufficient exposed surface')
    primary=min(fits,key=lambda f:abs(f['cutoff_dx']-1.5))
    return dict(status='estimated' if primary['optimizer_success'] and not primary['at_parameter_bound'] and primary['fit_rms']<.5*dx else 'unresolved fit',
                **primary,cutoff_sensitivity=fits,
                angle_cutoff_range_deg=max(f['angle_deg'] for f in fits)-min(f['angle_deg'] for f in fits),
                limitation=('Single convex zero-gravity cap; interface offset fitted' if fit_shift
                            else 'Single convex zero-gravity cap; intersection constrained to the physical wall plane'))


def write_vtp(path,vertices,faces):
    with Path(path).open('w') as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData>')
        f.write(f'<Piece NumberOfPoints="{len(vertices)}" NumberOfPolys="{len(faces)}"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        np.savetxt(f,vertices,fmt='%.12g')
        f.write('</DataArray></Points><Polys><DataArray type="Int64" Name="connectivity" format="ascii">\n')
        np.savetxt(f,faces,fmt='%d')
        f.write('</DataArray><DataArray type="Int64" Name="offsets" format="ascii">\n')
        np.savetxt(f,np.arange(1,len(faces)+1)*3,fmt='%d')
        f.write('</DataArray></Polys></Piece></PolyData></VTKFile>\n')


def analyze_file(path,meta,grid_dx=.25,mesh=None):
    dx=meta['dx'];rho0=meta.get('rho0',1.)
    vtp=Path(path).suffix.lower()=='.vtp'
    reader=read_vtp_samples if vtp else read_samples
    xyz,rho,masses=reader(path,rho0*dx**3)
    vertices,faces,v_iso=reconstruct(xyz,masses,4*dx,grid_dx*dx,rho0)
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    adjacency=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(vertices),len(vertices)))
    components,_=connected_components(adjacency,directed=False)
    volume=float(np.sum(masses/rho))
    mismatch=v_iso/volume-1
    result=dict(method='Wendland liquid density rho0/2 isosurface; volume-constrained zero-gravity cap',
                grid_dx=grid_dx,liquid_volume=volume,isosurface_volume=v_iso,
                isosurface_volume_relative_error=mismatch,volume_criterion_passed=abs(mismatch)<=.1,
                surface_vertices=len(vertices),surface_components=int(components),equilibrium_established=False)
    if meta['case']=='sessile':
        center=np.average(xyz[:,:2],axis=0,weights=masses)
        wall_z=meta.get('wall_z',0.)
        if vtp:
            # Scale the optimization to dx so the stopping tolerance also works in SI.
            origin=np.array([*center,wall_z])
            fit=fit_cap((vertices-origin)/dx,volume/dx**3,1.,np.zeros(2),fit_shift=False)
            for row in [fit,*fit.get('cutoff_sensitivity',[])]:
                for key in ('interface_shift','radius','center_z','fit_rms'):
                    if key in row:row[key]*=dx
                if 'center_z' in row:row['center_z']+=wall_z
            result['cap_fit']=fit
            result['touches_wall']=bool(xyz[:,2].min()-wall_z<=dx)
        else:
            result['cap_fit']=fit_cap(vertices,volume,dx,center)
        result['calibration_eligible']=(components==1 and abs(mismatch)<=.1 and
            result['cap_fit']['status']=='estimated' and result.get('touches_wall',True))
    if mesh: write_vtp(mesh,vertices,faces)
    return result


def history_summary(history):
    """Describe observed late-time drift; never declare equilibrium from a fit alone."""
    usable=[r for r in history if r.get('calibration_eligible') and
            r.get('angle_deg') is not None and math.isfinite(float(r['angle_deg']))]
    result=dict(frames=len(history),usable_frames=len(usable),equilibrium_established=False)
    if not history: return result
    end=max(float(r['time']) for r in history)
    tail=[r for r in usable if float(r['time'])>=.75*end]
    result['tail_usable_frames']=len(tail)
    if len(tail)>=3:
        t=np.array([float(r['time']) for r in tail])
        angles=np.array([float(r['angle_deg']) for r in tail])
        if np.ptp(t)>0:
            result.update(tail_start=float(t.min()),tail_end=float(t.max()),
                tail_angle_range_deg=float(np.ptp(angles)),
                tail_angle_slope_deg_per_time=float(np.polyfit(t-t.min(),angles,1)[0]),
                last_angle_deg=float(angles[-1]))
    return result


def analyze_folder(folder,grid_dx=.25,mesh=False):
    folder=Path(folder);meta=json.loads((folder/'pf_run.json').read_text())
    result=analyze_file(folder/'pf_final_liquid.csv',meta,grid_dx,folder/'pf_interface.vtp' if mesh else None)
    index=folder/'pf_snapshots.csv'
    if index.exists():
        with index.open() as f: frames=list(csv.DictReader(f))
        history=[]
        for frame in frames:
            data=result if int(frame['step'])==meta.get('steps') else analyze_file(folder/frame['file'],meta,grid_dx)
            fit=data.get('cap_fit',{})
            history.append(dict(step=frame['step'],time=frame['time'],angle_deg=fit.get('angle_deg'),
                        interface_shift=fit.get('interface_shift'),fit_rms=fit.get('fit_rms'),
                        volume_error=data['isosurface_volume_relative_error'],calibration_eligible=data.get('calibration_eligible',False)))
        if history:
            with (folder/'pf_angle_history.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(history[0]));writer.writeheader();writer.writerows(history)
        result['angle_history']=history_summary(history)
    (folder/'pf_interface_analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--grid-dx',type=float,default=.25,help='Grid spacing / sample spacing; paper uses 0.1')
    p.add_argument('--mesh',action='store_true',help='Write final fluid density isosurface as VTP')
    a=p.parse_args()
    if not math.isfinite(a.grid_dx) or not .05<=a.grid_dx<=1: p.error('Use 0.05 <= grid-dx <= 1')
    print(json.dumps(analyze_folder(a.directory,a.grid_dx,a.mesh),indent=2))
