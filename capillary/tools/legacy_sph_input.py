"""Normalize the original soot-dem/LB XML layout for the current SPH runner.

Preserve SI geometry and DEM coefficients. Never infer a calibrated sff/sfs from
an LB contact angle or target surface tension. No force model is changed here.
"""
import math
from pathlib import Path
import xml.etree.ElementTree as ET

UNITS={'m':1.,'mm':1e-3,'um':1e-6,'nm':1e-9}


def attrs(e,allowed):
    unknown=set(e.attrib)-set(allowed.split())
    if unknown:raise ValueError(f'Unknown {e.tag} attributes: {sorted(unknown)}')


def unit(value):
    if value not in UNITS:raise ValueError(f'Unsupported length_unit {value!r}; use m/mm/um/nm')
    return UNITS[value]


def real(value,name,positive=False):
    try:x=float(value)
    except (ValueError,TypeError) as exc:raise ValueError(f'Invalid {name}: {value}') from exc
    if not math.isfinite(x) or (positive and x<=0):raise ValueError(f'Invalid {name}: {value}')
    return x


def integer(value,name,minimum=0):
    try:x=int(value)
    except (ValueError,TypeError) as exc:raise ValueError(f'Invalid {name}: {value}') from exc
    if str(x)!=str(value).strip() or x<minimum:raise ValueError(f'Invalid integer {name}: {value}')
    return x


def declarations(root,path):
    values={};sources=[]
    def collect(element,where,stack):
        where=Path(where).resolve()
        if where in stack:raise ValueError(f'Cyclic XML include: {where}')
        if len(stack)>32:raise ValueError('XML include nesting exceeds 32')
        sources.append(str(where))
        for item in element.findall('let'):
            attrs(item,'id type');key=item.attrib['id'];kind=item.attrib['type'];text=(item.text or '').strip()
            if not key or key in values:raise ValueError(f'Duplicate/empty let id: {key}')
            if kind=='real':value=real(text,key)
            elif kind=='integer':
                try:value=int(text)
                except ValueError as exc:raise ValueError(f'Invalid integer {key}') from exc
            elif kind=='path':value=str((where.parent/text).resolve())
            elif kind=='string':value=text
            else:raise ValueError(f'Unsupported let type: {kind}')
            values[key]=(kind,value)
        for item in element.findall('include'):
            attrs(item,'');target=(where.parent/(item.text or '').strip()).resolve()
            child=ET.parse(target).getroot()
            if child.tag!='include_file':raise ValueError(f'Include root must be include_file: {target}')
            if any(x.tag not in ('let','include') for x in child):raise ValueError(f'Unexpected include-file element: {target}')
            collect(child,target,stack+[where])
    collect(root,path,[])
    return values,sources


def normalize(path,root,soot_defaults):
    attrs(root,'type')
    if root.get('type') not in ('capillary_dem','capillary_sph_dem'):
        raise ValueError('Simulation type must be capillary_sph_dem (or legacy capillary_dem)')
    allowed={'let','include','initialization','sph','domain_margin','output','lattice'}
    if any(e.tag not in allowed for e in root):raise ValueError('Unknown simulation XML section')
    for name in ('initialization','sph','domain_margin','output','lattice'):
        if len(root.findall(name))>1:raise ValueError(f'Duplicate {name} section')
    values,sources=declarations(root,path);used=set();notes=[]
    def get(name,default=None,required=False):
        if name not in values:
            if required:raise ValueError(f'Missing let parameter {name}')
            return default
        used.add(name);kind,value=values[name]
        if kind not in ('real','integer'):raise ValueError(f'{name} must be numeric')
        return value
    radius=real(get('r_part',required=True),'r_part',True)
    density=real(get('rho',required=True),'rho',True)
    dt=real(get('dt',required=True),'dt',True);end=real(get('t_tot',required=True),'t_tot',True)
    sph=root.find('sph')
    if sph is None:raise ValueError('Add an explicit <sph sound_speed="..." sff="..." sfs="..."> section; LB wetting parameters are not an SPH calibration')
    attrs(sph,'particle_spacing particles_per_primary_diameter length_unit sound_speed sff sfs dt_scale max_steps max_samples threads density_check fluid_equilibration_time ke_target_fraction ke_hold_time')
    for key in ('sound_speed','sff','sfs'):
        if key not in sph.attrib:raise ValueError(f'Missing explicit SPH parameter {key}')
    eq_let=get('fluid_equilibration_time')
    if eq_let is not None and sph.get('fluid_equilibration_time') is not None:raise ValueError('Specify fluid_equilibration_time once')
    eq_time=sph.get('fluid_equilibration_time',str(eq_let if eq_let is not None else 0))
    scale=unit(sph.get('length_unit','m'))
    cells=sph.get('particles_per_primary_diameter')
    old_cells=get('lb_cells_per_particle_diameter')
    if sph.get('particle_spacing') is not None:
        if cells is not None or old_cells is not None:raise ValueError('Choose particle_spacing or particles_per_primary_diameter, not both')
        dx=real(sph.get('particle_spacing'),'particle_spacing',True)*scale
    else:
        if cells is not None and old_cells is not None:raise ValueError('Specify resolution once, using the SPH field')
        cells=cells if cells is not None else old_cells
        if cells is None:raise ValueError('Specify SPH particle_spacing or particles_per_primary_diameter')
        if integer(str(cells),'particles_per_primary_diameter',1)!=float(cells):raise ValueError('Resolution must be integral')
        dx=2*radius/float(cells)
        if old_cells is not None:notes.append('Legacy LB cell count interpreted as SPH samples per primary diameter.')
    result=ET.Element('pf_dem')
    ET.SubElement(result,'fluid',dx=str(dx),rho0=str(real(get('liquid_density',required=True),'liquid_density',True)),
        sound_speed=str(real(sph.attrib['sound_speed'],'sound_speed',True)),
        viscosity=str(real(get('liquid_dynamic_viscosity',required=True),'liquid_dynamic_viscosity')),
        sff=str(real(sph.attrib['sff'],'sff')),sfs=str(real(sph.attrib['sfs'],'sfs')))
    init=root.find('initialization')
    if init is None:raise ValueError('Missing initialization section')
    attrs(init,'')
    if any(e.tag not in ('aggregate','liquid','necks') for e in init):raise ValueError('Unknown initialization element')
    for tag in ('aggregate','liquid','necks'):
        if len(init.findall(tag))>1:raise ValueError(f'Duplicate initialization {tag}')
    aggregate=init.find('aggregate')
    if aggregate is None:raise ValueError('Missing initialization aggregate')
    attrs(aggregate,'file format length_unit')
    if aggregate.get('format','vtk') not in ('vtk','vtp','csv'):raise ValueError('Aggregate format must be vtk/vtp/csv')
    # Original soot-dem VTK coordinates are in primary-radius units.
    raw_unit=aggregate.get('length_unit','primary_radius')
    aggregate_scale=radius if raw_unit=='primary_radius' else unit(raw_unit)
    ET.SubElement(result,'particles',file=str((path.parent/aggregate.attrib['file']).resolve()),
        radius=str(radius/aggregate_scale),density=str(density),length_scale=str(aggregate_scale),ignore_velocities='true')
    liquid=init.find('liquid')
    if liquid is None:raise ValueError('Missing initialization liquid')
    attrs(liquid,'mode alpha threshold rm length_unit')
    liquid_scale=unit(liquid.get('length_unit','m'))
    lm=liquid.get('mode','droplets')
    if lm=='levelset':
        if list(liquid):raise ValueError('Levelset liquid has no droplet children')
        ET.SubElement(result,'liquid',mode=lm,alpha=liquid.get('alpha','.5'),threshold=liquid.get('threshold','.5'),
                      rm=str(real(liquid.get('rm','0'),'rm')*liquid_scale))
    elif lm=='droplets':
        out=ET.SubElement(result,'liquid',mode=lm)
        for d in liquid:
            if d.tag!='droplet':raise ValueError('Expected droplet element')
            attrs(d,'x y z radius length_unit');ds=unit(d.get('length_unit',liquid.get('length_unit','m')))
            ET.SubElement(out,'droplet',**{k:str(real(d.get(k,'0') if k!='radius' else d.attrib[k],k,k=='radius')*ds) for k in ('x','y','z','radius')})
    else:raise ValueError('Liquid mode must be droplets or levelset')
    settings={}
    optional={'normal_yield':0,'shear_yield':0}
    for name,default in soot_defaults:
        if name=='max_dt':settings[name]=str(dt)
        elif name in optional:settings[name]=str(get(name,optional[name]))
        else:settings[name]=str(get(name,required=True))
    for name in ('neck_radius','normal_strength','shear_strength'):
        value=get(name)
        if value is not None:
            settings[name]=str(value)
    ET.SubElement(result,'soot_dem',**settings)
    neck=init.find('necks')
    if neck is None:ET.SubElement(result,'necks',mode='auto')
    else:
        attrs(neck,'mode file');kw=dict(neck.attrib)
        if 'file' in kw:kw['file']=str((path.parent/kw['file']).resolve())
        ET.SubElement(result,'necks',**kw)
    margin=root.find('domain_margin');kwargs={}
    if margin is not None:
        attrs(margin,'x y z length_unit');ms=unit(margin.get('length_unit','m'))
        kwargs={f'margin_{k}':str(real(margin.get(k,'0'),'margin_'+k)*ms) for k in 'xyz'}
    ET.SubElement(result,'domain',**kwargs,max_samples=sph.get('max_samples','2000000'))
    output=root.find('output');settings={}
    fields=('surface_every_sph_steps','fluid_every_sph_steps','particle_every_sph_steps','dem_every_sph_steps','log_every_sph_steps','checkpoint_every_sph_steps',
            'surface_every_lb_steps','fluid_every_lb_steps','particle_every_lb_steps','dem_every_dem_steps','log_every_lb_steps')
    if output is not None:
        attrs(output,' '.join(fields));settings=dict(output.attrib)
    def cadence(new,old,default):
        if new in settings and old in settings:raise ValueError(f'Specify {new} once, without its LB alias')
        value=settings.get(new,settings.get(old,str(default)))
        return integer(value,new)
    fluid_dump=cadence('fluid_every_sph_steps','fluid_every_lb_steps',100)
    particle_dump=cadence('particle_every_sph_steps','particle_every_lb_steps',fluid_dump)
    surface_dump=cadence('surface_every_sph_steps','surface_every_lb_steps',fluid_dump)
    record=cadence('log_every_sph_steps','log_every_lb_steps',100)
    dem_record=cadence('dem_every_sph_steps','dem_every_dem_steps',record)
    if fluid_dump!=particle_dump or dem_record!=record:
        raise ValueError('Current driver uses one fluid/particle dump interval and one DEM/log interval; set each pair equal')
    if surface_dump not in (0,fluid_dump):raise ValueError('Offline surface cadence must equal the saved fluid cadence (or zero)')
    if record==0:raise ValueError('Log/DEM record interval must be positive')
    run=ET.SubElement(result,'run',end_time=str(end),dt_scale=sph.get('dt_scale','.75'),
        max_steps=sph.get('max_steps','10000000'),record_every_steps=str(record),dump_every_steps=str(fluid_dump),
        checkpoint_every_steps=str(integer(settings.get('checkpoint_every_sph_steps','2000'),'checkpoint interval')),
        threads=sph.get('threads','1'),density_check=sph.get('density_check','true'),fluid_equilibration_time=eq_time)
    for key in ('ke_target_fraction','ke_hold_time'):
        value=get(key)
        if value is not None and sph.get(key) is not None:
            raise ValueError(f'Specify {key} once')
        if sph.get(key) is not None:
            run.set(key,sph.attrib[key])
        elif value is not None:
            run.set(key,str(value))
    targets={}
    for key in ('liquid_surface_tension','liquid_contact_angle_deg'):
        value=get(key)
        if value is not None:targets[key]=value
    for name in ('r_verlet','neighbor_update_period'):
        if name in values:used.add(name);notes.append(f'{name} is obsolete: cached SPH neighbors and native all-pair DEM forces do not use the old LB/DEM neighbor setting.')
    lattice=root.find('lattice')
    if lattice is not None:
        attrs(lattice,'initialization_samples_per_axis')
        integer(lattice.get('initialization_samples_per_axis','1'),'initialization_samples_per_axis',1)
        notes.append('Legacy VOF initialization samples are not used for SPH particle generation.')
    unknown=set(values)-used
    if unknown:raise ValueError(f'Unsupported let parameters: {sorted(unknown)}')
    notes.append('Surface meshes are exported offline from saved liquid frames; no online surface meshing is enabled.')
    return result,dict(format='soot-dem simulation XML',includes=sources,target_physics=targets,
                       calibration_status='provisional_explicit_sff_sfs',surface_every_sph_steps=surface_dump,notes=notes)
