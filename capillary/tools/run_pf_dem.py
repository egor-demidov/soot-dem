#!/usr/bin/env python3
"""Prepare XML aggregate/liquid geometry and run reciprocal DEM/PF-SPH coupling.

No calibration is inferred. Existing output directories are never overwritten.
Particle-file velocities are ignored by default. All input uses one consistent
unit system; this runner does not convert provisional PF strengths into SI gamma.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


def number(value, name, positive=False, nonnegative=False):
    try:x=float(value)
    except (TypeError,ValueError) as exc:raise ValueError(f'Invalid {name}: {value}') from exc
    if not math.isfinite(x) or (positive and x<=0) or (nonnegative and x<0):
        raise ValueError(f'Invalid {name}: {value}')
    return x


def integer(value, name, minimum=0):
    x=int(value)
    if str(x)!=str(value).strip() or x<minimum:
        raise ValueError(f'Invalid integer {name}: {value}')
    return x


def boolean(value):
    if value not in ('true','false','1','0'):
        raise ValueError(f'Expected true/false, got {value}')
    return value in ('true','1')


def check(element, allowed):
    unknown=set(element.attrib)-set(allowed.split())
    if unknown:raise ValueError(f'Unknown {element.tag} attributes: {sorted(unknown)}')


def vector(values, prefix='', default=0.):
    return [number(values.get(prefix+k,default),prefix+k) for k in 'xyz']


def particle(row, radius, density, scale=1., ignore_velocity=False):
    r=number(row.get('radius',radius),'radius',positive=True)*scale
    rho=number(row.get('density',density),'density',positive=True)
    fixed=boolean(str(row.get('fixed','false')).lower())
    velocity=[0.,0.,0.] if ignore_velocity else vector(row,'v')
    omega=[0.,0.,0.] if ignore_velocity else vector(row,'w')
    if fixed and any(velocity+omega):raise ValueError('Fixed particle has velocity/spin')
    return dict(center=[v*scale for v in vector(row)],radius=r,
                mass=4*math.pi*r**3*rho/3,fixed=fixed,velocity=velocity,omega=omega)


def vtk_rows(path, radius, density):
    if path.suffix.lower()=='.vtp':
        root=ET.parse(path).getroot();pieces=root.findall('.//Piece')
        if len(pieces)!=1:raise ValueError('Require a single-piece ASCII VTP aggregate')
        data=pieces[0].find('./Points/DataArray')
        if data is None or data.get('format','ascii')!='ascii':raise ValueError('Aggregate VTP Points must be ASCII')
        xyz=list(map(float,(data.text or '').split()))
        n=int(pieces[0].get('NumberOfPoints','0'))
        radii=None
        for a in pieces[0].findall('./PointData/DataArray'):
            if a.get('Name','').lower() in ('radius','radii'):
                if a.get('format','ascii')!='ascii':raise ValueError('Aggregate radius array must be ASCII')
                radii=list(map(float,(a.text or '').split()))
    else:
        text=path.read_text();lines=text.splitlines()
        if len(lines)<4 or lines[2].strip().upper()!='ASCII':
            raise ValueError('Aggregate VTK must be legacy ASCII; export ASCII from ParaView for this loader')
        tokens=' '.join(lines[3:]).split()
        try:
            i=tokens.index('POINTS');n=int(tokens[i+1]);start=i+3
            xyz=list(map(float,tokens[start:start+3*n]));radii=None
            if 'POINT_DATA' in tokens[start+3*n:]:
                p=tokens.index('POINT_DATA',start+3*n)
                if int(tokens[p+1])!=n:raise ValueError('VTK POINT_DATA count mismatch')
                k=p+2
                while k<len(tokens):
                    if tokens[k]=='SCALARS' and tokens[k+1].lower() in ('radius','radii'):
                        lookup=tokens.index('LOOKUP_TABLE',k);radii=list(map(float,tokens[lookup+2:lookup+2+n]));break
                    if tokens[k]=='CELL_DATA':break
                    k+=1
        except (IndexError,ValueError) as exc:
            raise ValueError('Malformed legacy VTK aggregate') from exc
    if n<1 or len(xyz)!=3*n or (radii is not None and len(radii)!=n):raise ValueError('Aggregate point/radius counts disagree')
    return [dict(x=xyz[3*i],y=xyz[3*i+1],z=xyz[3*i+2],
                 radius=radius if radii is None else radii[i],density=density) for i in range(n)]


SOOT_DEFAULTS=[('k_n',100),('gamma_n',.02),('k_t',20),('gamma_t',0),('mu_t',.5),('phi_t',.8),
 ('k_r',5),('gamma_r',0),('mu_r',.1),('phi_r',.8),('k_o',5),('gamma_o',0),('mu_o',.1),('phi_o',.8),
 ('k_n_bond',100),('gamma_n_bond',.02),('k_t_bond',20),('gamma_t_bond',0),
 ('k_r_bond',5),('gamma_r_bond',0),('k_o_bond',5),('gamma_o_bond',0),
 ('d_crit',.02),('A',0),('h0',.01),('max_dt',.0001),('normal_yield',0),('shear_yield',0)]

def neck_edges(path,n):
    path=Path(path)
    if path.suffix.lower()=='.csv':
        with path.open(newline='') as f:
            reader=csv.DictReader(f)
            if not {'i','j'}<=set(reader.fieldnames or []):raise ValueError('Neck CSV requires zero-based i,j columns')
            raw=[(integer(r['i'],'i'),integer(r['j'],'j')) for r in reader]
    elif path.suffix.lower()=='.vtk':
        lines=path.read_text().splitlines()
        if len(lines)<4 or lines[2].strip().upper()!='ASCII':raise ValueError('Neck VTK must be ASCII soot-dem neck output')
        tokens=' '.join(lines[3:]).split()
        try:
            k=tokens.index('connections')
            if int(tokens[k+1])!=2:raise ValueError('Neck connections must have two components')
            count=int(tokens[k+2])*2
            ids=[int(v) for v in tokens[k+4:k+4+count]]
            if len(ids)!=count or len(set(ids))!=len(ids):raise ValueError('Truncated/duplicate neck connection array')
            directed=set((q//n,q%n) for q in ids)
            if any((j,i) not in directed for i,j in directed):raise ValueError('Neck graph must be symmetric')
            raw=[(i,j) for i,j in directed if i<j]
            if any(i==j or i>=n or j>=n for i,j in directed):raise ValueError('Invalid neck index')
        except (ValueError,IndexError) as exc:raise ValueError('Invalid soot-dem neck VTK connections') from exc
    else:raise ValueError('Neck file must be CSV or ASCII soot-dem neck VTK')
    edges=[tuple(sorted(pair)) for pair in raw]
    if len(set(edges))!=len(edges) or any(i==j or i<0 or j>=n for i,j in edges):raise ValueError('Duplicate/invalid neck indices')
    return sorted(edges)


def parse_xml(path):
    path=Path(path).resolve();root=ET.parse(path).getroot()
    compatibility=None
    if root.tag=='simulation':
        from legacy_sph_input import normalize
        root,compatibility=normalize(path,root,SOOT_DEFAULTS)
    if root.tag!='pf_dem':raise ValueError('Root element must be pf_dem or simulation')
    check(root,'')
    known={'fluid','contacts','particles','particle','liquid','domain','run','soot_dem','necks','wall'}
    for child in root:
        if child.tag not in known:raise ValueError(f'Unknown XML element {child.tag}')
        if child.tag!='liquid' and list(child):raise ValueError(f'Unexpected nested element in {child.tag}')
        if child.tag not in ('particle',) and len(root.findall(child.tag))!=1:
            raise ValueError(f'Duplicate XML section {child.tag}')
    f=root.find('fluid')
    if f is None:raise ValueError('Missing fluid section')
    check(f,'dx rho0 sound_speed viscosity sff sfs')
    fluid=[number(f.get(k,default),k,positive=k in ('dx','rho0','sound_speed'),nonnegative=k in ('viscosity','sff','sfs'))
           for k,default in [('dx','.08'),('rho0','1'),('sound_speed','10'),('viscosity','.05'),('sff','.001'),('sfs','.003')]]
    c=root.find('contacts')
    if c is not None:check(c,'stiffness damping friction slip_speed')
    contact=[number(c.get(k,d) if c is not None else d,k,positive=k=='slip_speed',nonnegative=k!='slip_speed')
             for k,d in [('stiffness','100'),('damping','0'),('friction','0'),('slip_speed','.01')]]
    bodies=[];aggregate=root.find('particles')
    if aggregate is not None:
        check(aggregate,'file radius density length_scale ignore_velocities')
        source=(path.parent/aggregate.attrib['file']).resolve()
        radius=aggregate.get('radius');density=aggregate.get('density','2')
        scale=number(aggregate.get('length_scale','1'),'length_scale',positive=True)
        if source.suffix.lower() in ('.vtk','.vtp'):rows=vtk_rows(source,radius,density)
        elif source.suffix.lower()=='.csv':
            with source.open(newline='') as file:
                reader=csv.DictReader(file)
                if not {'x','y','z'}<=set(reader.fieldnames or []):raise ValueError('Aggregate CSV requires x,y,z')
                rows=list(reader)
        else:raise ValueError('Particle file must be .csv, legacy ASCII .vtk, or ASCII .vtp')
        for row in rows:
            if radius is None and not row.get('radius'):raise ValueError('Specify a radius column/array or particles radius')
            bodies.append(particle(row,radius,density,scale,boolean(aggregate.get('ignore_velocities','true'))))
    for e in root.findall('particle'):
        check(e,'x y z radius density fixed vx vy vz wx wy wz')
        bodies.append(particle(e.attrib,e.get('radius','.24'),e.get('density','2')))
    l=root.find('liquid')
    if l is None:raise ValueError('Missing liquid section')
    check(l,'mode alpha threshold rm')
    mode=l.get('mode','droplets')
    if mode not in ('droplets','levelset'):raise ValueError('Liquid mode must be droplets or levelset')
    alpha=number(l.get('alpha','.5'),'alpha',positive=True)
    threshold=number(l.get('threshold','.5'),'threshold',positive=True)
    rm=number(l.get('rm','0'),'rm',nonnegative=True)
    drops=[]
    for e in l:
        if e.tag!='droplet':raise ValueError('Expected droplet element inside liquid')
        check(e,'x y z radius');drops.append(dict(center=vector(e.attrib),radius=number(e.attrib['radius'],'droplet radius',positive=True)))
    if mode=='droplets' and not drops:raise ValueError('No liquid droplets')
    if mode=='levelset' and (drops or threshold>=len(bodies)):
        raise ValueError('Levelset mode needs no droplets and threshold < particle count')
    wall_element=root.find('wall')
    wall=None
    if wall_element is not None:
        check(wall_element,'z half_width')
        if bodies:raise ValueError('Flat wall is only supported for sessile-drop calibration')
        wall=dict(z=number(wall_element.get('z','0'),'wall z'),
                  half_width=number(wall_element.attrib['half_width'],'wall half_width',positive=True))
    if not bodies:
        if aggregate is not None or any(root.find(tag) is not None for tag in ('contacts','soot_dem','necks')):
            raise ValueError('Calibration inputs cannot specify DEM particles or contacts')
        if mode!='droplets' or len(drops)!=1:
            raise ValueError('Calibration requires exactly one droplet')
        if wall is not None:
            if abs(drops[0]['center'][2]-wall['z'])>=drops[0]['radius']:
                raise ValueError('Initial droplet must intersect the flat wall plane')
            if wall['half_width']<drops[0]['radius']+4*fluid[0]:
                raise ValueError('Wall half_width must cover the initial droplet radius plus 4 dx')
    domain=root.find('domain')
    if domain is not None:check(domain,'margin margin_x margin_y margin_z max_samples')
    margin=number(domain.get('margin',str(4*fluid[0])) if domain is not None else 4*fluid[0],'margin',nonnegative=True)
    max_samples=integer(domain.get('max_samples','2000000') if domain is not None else '2000000','max_samples',1)
    extents=[(b['center'],b['radius']) for b in bodies]+[(d['center'],d['radius']) for d in drops]
    if mode=='levelset':
        # If sum psi >= threshold then at least one psi >= threshold/N.
        # Bound its erfc tail by bisection, avoiding a SciPy dependency here.
        target=2*threshold/len(bodies);lo,hi=-30.,30.
        for _ in range(100):
            mid=.5*(lo+hi)
            if math.erfc(mid)>target:lo=mid
            else:hi=mid
        for b in bodies:extents.append((b['center'],max(b['radius'],(rm or b['radius'])*(1+alpha*hi))))
    margins=[number(domain.get('margin_'+axis,str(margin)) if domain is not None else margin,'margin_'+axis,nonnegative=True) for axis in 'xyz']
    lower=[min(x[k]-r for x,r in extents)-margins[k] for k in range(3)]
    upper=[max(x[k]+r for x,r in extents)+margins[k] for k in range(3)]
    run=root.find('run')
    if run is not None:check(run,'end_time dt_scale max_steps record_every_steps dump_every_steps checkpoint_every_steps threads density_check fluid_equilibration_time ke_target_fraction ke_hold_time')
    runtime={}
    for k,d in [('end_time','.5'),('dt_scale','.75')]:runtime[k]=number(run.get(k,d) if run is not None else d,k,positive=True)
    if runtime['dt_scale']>1:raise ValueError('dt_scale must be <=1')
    for k,d,minval in [('max_steps','1000000',1),('record_every_steps','100',1),('dump_every_steps','500',0),('checkpoint_every_steps','500',0),('threads','1',1)]:
        runtime[k]=integer(run.get(k,d) if run is not None else d,k,minval)
    runtime['fluid_equilibration_time']=number(run.get('fluid_equilibration_time','0') if run is not None else '0','fluid_equilibration_time',nonnegative=True)
    runtime['density_check']=boolean(run.get('density_check','true') if run is not None else 'true')
    if run is not None:
        if 'ke_target_fraction' in run.attrib:runtime['ke_target_fraction']=ke_fraction(run.attrib['ke_target_fraction'])
        if 'ke_hold_time' in run.attrib:runtime['ke_hold_time']=number(run.attrib['ke_hold_time'],'ke_hold_time',nonnegative=True)
    model=dict(fluid=fluid,contact=contact,bodies=bodies,mode=mode,alpha=alpha,threshold=threshold,rm=rm,drops=drops,lower=lower,upper=upper,max_samples=max_samples)
    if not bodies:
        if runtime['fluid_equilibration_time']:
            raise ValueError('fluid_equilibration_time applies to runs with DEM particles')
        model['calibration']='sessile' if wall is not None else 'droplet'
        model['wall']=wall
    soot=root.find('soot_dem');neck=root.find('necks')
    if soot is not None:
        strength_keys=('neck_radius','normal_strength','shear_strength')
        check(soot,' '.join([k for k,_ in SOOT_DEFAULTS]+list(strength_keys)))
        coefficients={k:number(soot.get(k,str(d)),k,positive=k in ('h0','max_dt'),nonnegative=True)
                      for k,d in SOOT_DEFAULTS}
        if any(key in soot.attrib for key in strength_keys):
            if not all(key in soot.attrib for key in strength_keys):
                raise ValueError('Specify neck_radius, normal_strength, and shear_strength together')
            if coefficients['normal_yield'] or coefficients['shear_yield']:
                raise ValueError('Neck strengths cannot be combined with nonzero normal_yield or shear_yield')
            radius,normal_strength,shear_strength=(
                number(soot.attrib[key],key,positive=True) for key in strength_keys)
            area=math.pi*radius*radius
            coefficients['normal_yield']=number(area*normal_strength,'normal yield force',positive=True)
            coefficients['shear_yield']=number(area*shear_strength,'shear yield force',positive=True)
        values=[coefficients[k] for k,_ in SOOT_DEFAULTS]
        if len(bodies)>2000:raise ValueError('soot-dem dense history arrays support at most 2000 particles in this driver')
        if c is not None and any(contact[:3]):raise ValueError('Do not combine built-in contacts with soot_dem')
        model['contact']=contact=[0.,0.,0.,.01]
        for b in bodies:
            for key in ('radius','mass'):
                if not math.isclose(b[key],bodies[0][key],rel_tol=1e-12,abs_tol=0):raise ValueError('soot-dem requires equal primary radii and masses')
        for k in (2,6,10):
            if values[k]==0 and (values[k+1]>0 or values[k+2]>0):raise ValueError('Positive shear stiffness required for friction/damping')
        mode_neck='auto';edges=[]
        if neck is not None:
            check(neck,'mode file')
            mode_neck=neck.get('mode','file' if neck.get('file') else 'auto')
            if mode_neck not in ('auto','none','file'):raise ValueError('Neck mode must be auto, none, or file')
            if mode_neck=='file':edges=neck_edges(path.parent/neck.attrib['file'],len(bodies))
            elif neck.get('file'):raise ValueError('Neck file requires mode=file')
        model['soot_dem']=dict(values=values,neck_mode={'none':0,'auto':1,'file':2}[mode_neck],edges=edges)
    elif neck is not None:raise ValueError('Necks require soot_dem force model')

    signature=hashlib.sha256(json.dumps(model,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    result=dict(model=model,run=runtime,signature=signature)
    if compatibility is not None:result['input_compatibility']=compatibility
    return result


def prepared_text(config):
    m=config['model'];r=config['run']
    rows=[['PFDEM_INPUT',3 if 'calibration' in m else (2 if 'soot_dem' in m else 1)],[config['signature']],m['fluid'],m['contact'],
          [r[k] for k in ('end_time','dt_scale','max_steps','record_every_steps','dump_every_steps','checkpoint_every_steps','threads')],
          m['lower']+m['upper'],[m['max_samples']],[len(m['bodies'])]]
    for b in m['bodies']:rows.append(b['center']+[b['radius'],b['mass'],int(b['fixed'])]+b['velocity']+b['omega'])
    rows += [[int(m['mode']=='levelset'),m['alpha'],m['threshold'],m['rm']],[len(m['drops'])]]
    for d in m['drops']:rows.append(d['center']+[d['radius']])
    if 'calibration' in m:
        wall=m['wall']
        rows.append([1,wall['z'],wall['half_width']] if wall is not None else [0])
    if 'soot_dem' in m:
        soot=m['soot_dem'];rows.extend([soot['values'],[soot['neck_mode'],len(soot['edges'])]])
        rows.extend(soot['edges'])
    return '\n'.join(' '.join(str(x) for x in row) for row in rows)+'\n'


def ke_fraction(value):
    fraction=number(value,'ke_target_fraction',positive=True)
    if fraction>=1:
        raise ValueError('ke_target_fraction must be strictly between 0 and 1')
    return fraction


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('xml',type=Path);p.add_argument('--executable',default='build/pf_dem')
    p.add_argument('--output',required=True,type=Path);p.add_argument('--resume',type=Path)
    p.add_argument('--end-time',type=float);p.add_argument('--steps',type=int)
    p.add_argument('--ke-target-fraction',type=float,help='Stop when coupled kinetic energy falls to this fraction of its running peak')
    p.add_argument('--ke-hold-time',type=float,metavar='SECONDS',help='Require kinetic energy to stay below the target for this duration')
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--disable-density-check',action='store_true',help='Disable only the reference-density stopping screen')
    p.add_argument('--legacy-overlap-mask',action='store_true',help='Preserve the old boundary rule when resuming an old overlapping checkpoint')
    p.add_argument('--add-liquid-volume',type=float,metavar='M3',help='Add a liquid dose on restart')
    p.add_argument('--liquid-addition-mode',choices=['liquid-vapor','uniform-liquid-vapor','uniform-solid-vapor','curvature-based'],
                   help='Placement rule for a dose (default: liquid-vapor, closest contour sites first)')
    p.add_argument('--saturation-ratio',type=float,help='Dimensionless saturation ratio for curvature-based addition')
    p.add_argument('--kelvin-length',type=float,metavar='METERS',help='Kelvin length for curvature-based addition')
    p.add_argument('--liquid-addition-seed',type=int,help='Random seed for curvature-based addition (default: 0)')
    p.add_argument('--post-addition-equilibration-time',type=float,metavar='SECONDS',help='Freeze DEM bodies for this duration after the dose')
    a=p.parse_args()
    try:
        config=parse_xml(a.xml)
        runtime=config['run']
        if a.ke_target_fraction is not None:runtime['ke_target_fraction']=ke_fraction(a.ke_target_fraction)
        if a.ke_hold_time is not None:runtime['ke_hold_time']=number(a.ke_hold_time,'ke_hold_time',nonnegative=True)
        if 'ke_hold_time' in runtime and 'ke_target_fraction' not in runtime:
            raise ValueError('ke_hold_time requires ke_target_fraction')
        if a.add_liquid_volume is not None:
            number(a.add_liquid_volume,'add-liquid-volume',positive=True)
            if a.resume is None:
                raise ValueError('--add-liquid-volume requires --resume')
            if not config['model']['bodies']:
                raise ValueError('Liquid addition requires DEM bodies')
            config['run']['add_liquid_volume']=a.add_liquid_volume
            config['run']['liquid_addition_mode']=a.liquid_addition_mode or 'liquid-vapor'
        elif a.liquid_addition_mode is not None:
            raise ValueError('--liquid-addition-mode requires --add-liquid-volume')
        curvature_based=a.liquid_addition_mode=='curvature-based'
        if curvature_based:
            if a.saturation_ratio is None or a.kelvin_length is None:
                raise ValueError('Curvature-based addition requires --saturation-ratio and --kelvin-length')
            number(a.saturation_ratio,'saturation-ratio',positive=True)
            number(a.kelvin_length,'kelvin-length',nonnegative=True)
            seed=0 if a.liquid_addition_seed is None else a.liquid_addition_seed
            if not 0<=seed<2**64:
                raise ValueError('Liquid addition seed must be an unsigned 64-bit integer')
            config['run'].update(saturation_ratio=a.saturation_ratio,kelvin_length_m=a.kelvin_length,
                                 liquid_addition_seed=seed)
        elif any(value is not None for value in (a.saturation_ratio,a.kelvin_length,a.liquid_addition_seed)):
            raise ValueError('Kelvin inputs and the addition seed require --liquid-addition-mode curvature-based')
        if a.post_addition_equilibration_time is not None:
            number(a.post_addition_equilibration_time,'post-addition-equilibration-time',nonnegative=True)
            if a.add_liquid_volume is None:
                raise ValueError('--post-addition-equilibration-time requires --add-liquid-volume')
            config['run']['post_addition_equilibration_time']=a.post_addition_equilibration_time
        if a.disable_density_check:config['run']['density_check']=False
        if a.output.exists():raise ValueError('Output exists; choose a fresh directory (including for restart)')
        if a.resume and not a.resume.is_file():raise ValueError('Checkpoint not found')
        exe=shutil.which(a.executable) or (str(Path(a.executable).resolve()) if Path(a.executable).is_file() else None)
        if exe:exe=str(Path(exe).resolve())
        if not exe and not a.prepare_only:raise ValueError('Build pf_dem or provide --executable')
        if a.end_time is not None:number(a.end_time,'end-time',positive=True)
        if a.steps is not None and a.steps<1:raise ValueError('steps must be positive')
        a.output.mkdir(parents=True)
        prepared=(a.output/'prepared.dat').resolve();prepared.write_text(prepared_text(config))
        (a.output/'resolved_input.json').write_text(json.dumps(config,indent=2)+'\n')
        shutil.copyfile(a.xml,a.output/'input.xml')
        cmd=[exe or a.executable,'--input',str(prepared)]
        if not config['run'].get('density_check',True):cmd+=['--disable-density-check']
        if a.resume is None and config['run'].get('fluid_equilibration_time',0)>0:cmd+=['--fluid-equilibration-time',str(config['run']['fluid_equilibration_time'])]
        if 'ke_target_fraction' in runtime:cmd+=['--ke-target-fraction',str(runtime['ke_target_fraction'])]
        if 'ke_hold_time' in runtime:cmd+=['--ke-hold-time',str(runtime['ke_hold_time'])]
        if a.legacy_overlap_mask:cmd+=['--legacy-overlap-mask']
        if a.resume:cmd+=['--resume',str(a.resume.resolve())]
        if a.add_liquid_volume is not None:cmd+=['--add-liquid-volume',str(a.add_liquid_volume)]
        if a.liquid_addition_mode is not None:cmd+=['--liquid-addition-mode',a.liquid_addition_mode]
        if curvature_based:
            cmd+=['--saturation-ratio',str(a.saturation_ratio),'--kelvin-length',str(a.kelvin_length),
                  '--liquid-addition-seed',str(seed)]
        if a.post_addition_equilibration_time is not None:
            cmd+=['--post-addition-equilibration-time',str(a.post_addition_equilibration_time)]
        if a.end_time is not None:cmd+=['--end-time',str(a.end_time)]
        if a.steps is not None:cmd+=['--steps',str(a.steps)]
        (a.output/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        print(shlex.join(cmd),flush=True)
        if a.prepare_only:return 0
        with (a.output/'run.log').open('w') as log:
            result=subprocess.run(cmd,cwd=a.output,stdout=log,stderr=subprocess.STDOUT)
        print(f'Exit {result.returncode}; output: {a.output}; log: {a.output/"run.log"}')
        if result.returncode:print((a.output/'run.log').read_text()[-4000:],file=sys.stderr)
        return result.returncode
    except (ValueError,KeyError,OSError,ET.ParseError) as exc:p.error(str(exc))

if __name__=='__main__':raise SystemExit(main())
