#!/usr/bin/env python3
"""Export DEM liquid VTP point clouds as Wendland-density surface meshes.

Postprocessing only: no new dynamics, density resets, volume fitting or smoothing
of the simulated coordinates. Use opaque particle meshes to hide portions of the
reconstructed liquid-density surface inside the solid geometry.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
from pf_interface import density_grid


def read_liquid(path):
    root=ET.parse(path).getroot();pieces=root.findall('.//Piece')
    if len(pieces)!=1:raise ValueError(f'Expected one VTP piece: {path}')
    piece=pieces[0];a=piece.find('./Points/DataArray')
    if a is None or a.get('format','ascii')!='ascii':raise ValueError('Expected ASCII point coordinates from pf_dem')
    points=np.fromstring(a.text or '',sep=' ')
    count=int(piece.get('NumberOfPoints','-1'))
    if count<1 or points.size!=3*count or not np.isfinite(points).all():raise ValueError(f'Invalid liquid points: {path}')
    time=None
    for a in root.findall('.//FieldData/DataArray'):
        if a.get('Name')=='TimeValue' and a.get('format','ascii')=='ascii':
            time=float(a.text)
            if not math.isfinite(time):raise ValueError('Nonfinite frame time')
    return points.reshape((-1,3)),time


def write_surface(path,vertices,faces,normals):
    with Path(path).open('w') as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData>')
        f.write(f'<Piece NumberOfPoints="{len(vertices)}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{len(faces)}">')
        f.write('<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        np.savetxt(f,vertices,fmt='%.17g')
        f.write('</DataArray></Points><PointData Normals="Normals"><DataArray type="Float64" Name="Normals" NumberOfComponents="3" format="ascii">\n')
        np.savetxt(f,normals,fmt='%.9g')
        f.write('</DataArray></PointData><Polys><DataArray type="Int64" Name="connectivity" format="ascii">\n')
        np.savetxt(f,faces,fmt='%d')
        f.write('</DataArray><DataArray type="Int64" Name="offsets" format="ascii">\n')
        np.savetxt(f,3*np.arange(1,len(faces)+1),fmt='%d')
        f.write('</DataArray></Polys></Piece></PolyData></VTKFile>\n')


def export(directory,output,grid_dx=.5,every=1,last_only=False,max_nodes=30000000):
    # Import before creating output so missing dependencies leave no partial folder.
    from skimage.measure import marching_cubes
    directory=Path(directory);output=Path(output)
    if output.exists():raise ValueError('Output directory exists; select a fresh --output')
    if not (grid_dx>0 and math.isfinite(grid_dx)) or every<1 or max_nodes<8:raise ValueError('Invalid grid spacing/frame stride/node budget')
    meta=json.loads((directory/'resolved_input.json').read_text())
    dx,rho0=map(float,meta['model']['fluid'][:2])
    if not all(math.isfinite(v) and v>0 for v in (dx,rho0)):raise ValueError('Invalid fluid dx/rho0')
    spacing=grid_dx*dx;support=4*dx;mass=rho0*dx**3
    files=[]
    for path in directory.glob('liquid_*.vtp'):
        match=re.fullmatch(r'liquid_(\d+)\.vtp',path.name)
        if match:files.append((int(match[1]),path))
    files.sort()
    if not files:raise ValueError('No liquid point-cloud frames; enable dump_every_steps in the simulation XML')
    selected=files[-1:] if last_only else files[::every]
    if not last_only and selected[-1]!=files[-1]:selected.append(files[-1])
    times={}
    history=directory/'dem_history.csv'
    if history.is_file():
        with history.open(newline='') as f:
            for row in csv.DictReader(f):times[int(row['step'])]=float(row['time'])
    output.mkdir(parents=True);reports=[]
    for step,path in selected:
        xyz,time=read_liquid(path)
        if time is None:time=times.get(step)
        grid,origin=density_grid(xyz,np.full(len(xyz),mass),support,spacing,max_nodes=max_nodes)
        level=rho0/2
        resolved=bool(grid.min()<level<grid.max())
        if resolved:
            v,faces,normals,_=marching_cubes(grid,level=level,spacing=(spacing,)*3,allow_degenerate=False)
            v+=origin
        else:
            v=np.empty((0,3));faces=np.empty((0,3),dtype=int);normals=np.empty((0,3))
        name=f'liquid_surface_{step:010d}.vtp';write_surface(output/name,v,faces,normals)
        reports.append(dict(step=step,time=time,file=name,points=len(xyz),vertices=len(v),triangles=len(faces),resolved=resolved))
        print(f'{path.name}: {len(v)} surface vertices, {len(faces)} triangles'+('' if resolved else ' (unresolved/empty isosurface)'),flush=True)
    physical_times=all(r['time'] is not None and math.isfinite(r['time']) for r in reports)
    collection=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian')
    group=ET.SubElement(collection,'Collection')
    for r in reports:ET.SubElement(group,'DataSet',timestep=str(r['time'] if physical_times else r['step']),group='',part='0',file=r['file'])
    ET.ElementTree(collection).write(output/'liquid_surface.pvd',encoding='utf-8',xml_declaration=True)
    result=dict(method='Liquid mass-weighted Wendland kernel, support 4 dx, rho0/2 isosurface',
                grid_spacing=spacing,grid_dx=grid_dx,animation_time='simulation_time' if physical_times else 'step',
                physics_validated=False,frames=reports,
                limitations=['Visualization only; isosurface volume is not conserved liquid mass.',
                             'No solid clipping: display opaque particles to hide the reconstruction within solids.',
                             'Thin films, small gaps and contact-line shapes depend on particle resolution and kernel support.'])
    (output/'surface_export.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not physical_times:print('Some frame times unavailable; PVD animation uses step numbers for every frame.')
    print(f'Open in ParaView: {output/"liquid_surface.pvd"}')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--grid-dx',type=float,default=.5,help='Grid spacing / SPH particle spacing; .5 is a visual-quality starting point')
    p.add_argument('--every',type=int,default=1,help='Export every N saved frames; always include the last frame')
    p.add_argument('--last-only',action='store_true');p.add_argument('--max-nodes',type=int,default=30000000)
    a=p.parse_args()
    try:export(a.directory,a.output or a.directory/'liquid_surfaces',a.grid_dx,a.every,a.last_only,a.max_nodes)
    except ImportError:p.error('Install reconstruction dependencies: python3 -m pip install -r requirements.txt')
    except (ValueError,OSError,KeyError,ET.ParseError) as exc:p.error(str(exc))

if __name__=='__main__':main()
