#!/usr/bin/env python3
"""Prepare a build-local soot-dem header overlay; never edit the user's checkout.

Adds history/timestep accessors required for transactional coupling and checkpoint
I/O. Corrects torsion-history projection onto the axial rather than tangent plane.
Tested soot-dem 88ae15b, libgran ad9eb4d. Fails on an unrecognized API.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def replace_once(text,old,new):
    if text.count(old)!=1:raise ValueError('Unsupported soot-dem/libgran API; expected header pattern not found exactly once: '+old[:70])
    return text.replace(old,new,1)


def prepare(source,output):
    source=Path(source).resolve();output=Path(output).resolve()
    for protected in (source,source/'src',source/'deps'):
        if output==protected or (protected!=source and protected in output.parents):
            raise ValueError('Overlay must not replace source directories')
    names=['src/aggregate.h','src/break_neck.h','deps/libgran/include/libgran/contact_force/contact_force.h',
           'deps/libgran/include/libgran/sinter_bridge/alt_sinter_bridge.h']
    for name in names:
        if not (source/name).is_file():raise ValueError(f'Missing {source/name}; initialize libgran and Eigen submodules')
    if not (source/'deps/eigen/Eigen/Core').is_file():raise ValueError('Missing soot-dem Eigen submodule')
    output.mkdir(parents=True,exist_ok=True)
    aggregate=(source/names[0]).read_text()
    aggregate=replace_once(aggregate,'#include <libgran/granular_system/granular_system.h>','#include <libgran/utility/operator_overloads.h>')
    aggregate=replace_once(aggregate,'private:', '''
    sinter_model_t & pf_sinter_model() { return sinter_model; }
    sinter_model_t const & pf_sinter_model() const { return sinter_model; }
private:''')
    (output/'aggregate.h').write_text(aggregate)
    shutil.copy2(source/'src/break_neck.h',output/'break_neck.h')
    for relative in names[2:]:
        p=source/relative;s=p.read_text()
        s=replace_once(s,'r_part, mass, inertia, dt, real_zero;', 'r_part, mass, inertia;\n    real_t dt;\n    const real_t real_zero;')
        # The torsion spring is axial; sliding/rolling springs are tangential.
        s=replace_once(s,'field_value_t xi_new = xi - xi.dot(n) * n;',
                      'field_value_t xi_new;\n        if constexpr (model_num == 2) xi_new = xi.dot(n) * n;\n        else xi_new = xi - xi.dot(n) * n;')
        methods='''
    void pf_set_dt(real_t value) { dt = value; PF_FORWARD_DT }
    auto & pf_springs() { return contact_springs; }
    auto const & pf_springs() const { return contact_springs; }
'''
        if 'alt_sinter_bridge' in relative:
            methods=methods.replace('PF_FORWARD_DT','contact_force.pf_set_dt(value);')
            methods+='''    auto & pf_nonbonded() { return contact_force; }
    auto const & pf_nonbonded() const { return contact_force; }
'''
        else:methods=methods.replace('PF_FORWARD_DT','')
        s=replace_once(s,'private:',methods+'\nprivate:')
        dest=output/Path(relative).relative_to('deps/libgran/include')
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(s)
    try:commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    except (OSError,subprocess.CalledProcessError):commit='unknown'
    tracked=[source/'src/aggregate.h',source/'src/break_neck.h']+sorted((source/'deps/libgran/include').rglob('*.h'))
    digest=hashlib.sha256()
    for file in tracked:
        digest.update(str(file.relative_to(source)).encode());digest.update(file.read_bytes())
    digest.update(Path(__file__).read_bytes())
    digest.update((Path(__file__).resolve().parents[1]/'src/pf_soot_dem.h').read_bytes())
    fingerprint=digest.hexdigest()
    (output/'pf_soot_build_id.h').write_text('#pragma once\n#define PF_SOOT_BACKEND_FINGERPRINT "'+fingerprint+'"\n')
    (output/'backend.json').write_text(json.dumps(dict(soot_dem_source=str(source),commit=commit,fingerprint=fingerprint,
       modifications=['transactional history access','variable spring update dt','axial torsion-history projection'],
       source_checkout_modified=False),indent=2)+'\n')
    print(output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    try:prepare(a.source,a.output)
    except (ValueError,OSError) as exc:p.error(str(exc))
