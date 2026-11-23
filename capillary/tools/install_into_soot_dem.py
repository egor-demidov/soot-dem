#!/usr/bin/env python3
"""Merge the capillary component into soot-dem without replacing its target list.

If the root CMakeLists.txt was overwritten by the standalone capillary build,
recover the most recent soot-dem version from its local Git history. Back up the
current root file first. Preserve original sources, dependencies, runs and builds.
"""
import argparse
import json
from pathlib import Path,PurePosixPath
import re
import shutil
import subprocess

INCLUDE='include("${CMAKE_CURRENT_LIST_DIR}/cmake/SootCapillaryCoupling.cmake")'


def soot_build(text):
    return bool(re.search(r'add_executable\s*\(\s*restructuring\b',text) and
                re.search(r'add_executable\s*\(\s*aggregation\b',text))


def recover_original(repo,current):
    if soot_build(current):return current,'working tree'
    if 'sphinxsys_capillary_validation' not in current:
        raise ValueError('Root CMake is not a recognized soot-dem or overwritten capillary build; refusing to replace it')
    def git(*args):
        return subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True,check=True).stdout
    try:
        # HEAD can retain a user's custom target definitions. If it was committed
        # after the accidental replacement, search older root-file revisions.
        head=git('show','HEAD:CMakeLists.txt')
        if soot_build(head):return head,'Git HEAD'
        commits=git('log','-n','100','--format=%H','--','CMakeLists.txt').splitlines()
        for commit in commits:
            text=git('show',commit+':CMakeLists.txt')
            if soot_build(text):return text,'Git '+commit
    except (OSError,subprocess.CalledProcessError) as exc:
        raise ValueError('Cannot recover the original soot-dem CMake from local Git; restore that file first, then rerun this installer') from exc
    raise ValueError('No original soot-dem CMake found in the last 100 file revisions; restore it first')


def install(repo,source=None):
    repo=Path(repo).resolve()
    source=Path(source).resolve() if source else Path(__file__).resolve().parents[1]
    cmake=repo/'CMakeLists.txt'
    for file in ('src/aggregate.h','src/restructuring.cpp','src/aggregation.cpp'):
        if not (repo/file).is_file():raise ValueError(f'Not a soot-dem source checkout: missing {repo/file}')
    if not cmake.is_file():raise ValueError('Original root CMakeLists.txt is missing')
    current=cmake.read_text();original,origin=recover_original(repo,current)
    manifest=json.loads((source/'cmake/capillary_payload.json').read_text())
    if not isinstance(manifest,list) or not manifest:raise ValueError('Invalid coupling payload manifest')
    # Resolve every input before changing the destination; paths are constrained
    # beneath the package and the dedicated capillary component.
    plan=[]
    for entry in manifest:
        src=PurePosixPath(entry['source']);dst=PurePosixPath(entry['destination'])
        if src.is_absolute() or dst.is_absolute() or '..' in src.parts or '..' in dst.parts:
            raise ValueError('Unsafe coupling payload path')
        actual=source/str(src)
        if not actual.is_file():raise ValueError(f'Incomplete coupling package: {actual}')
        plan.append((actual,repo/'capillary'/str(dst)))
    if not (source/'cmake/SootCapillaryCoupling.cmake').is_file():raise ValueError('Missing root integration module')
    updated=original.rstrip()+'\n'
    if INCLUDE not in updated:
        updated+='\n# PF-SPH fluid coupling, in addition to all original soot-dem targets.\n'+INCLUDE+'\n'
    backup=None
    if updated!=current:
        backup=repo/'CMakeLists.txt.before-capillary-merge'
        i=1
        while backup.exists():
            backup=repo/f'CMakeLists.txt.before-capillary-merge.{i}';i+=1
        shutil.copy2(cmake,backup)
    for src,dst in plan:
        dst.parent.mkdir(parents=True,exist_ok=True)
        if src.resolve()!=dst.resolve():shutil.copy2(src,dst)
    module=repo/'cmake/SootCapillaryCoupling.cmake';module.parent.mkdir(parents=True,exist_ok=True)
    if (source/'cmake/SootCapillaryCoupling.cmake').resolve()!=module.resolve():shutil.copy2(source/'cmake/SootCapillaryCoupling.cmake',module)
    if updated!=current:
        temporary=repo/'CMakeLists.txt.capillary-merge.tmp';temporary.write_text(updated);temporary.replace(cmake)
    result=dict(repository=str(repo),original_cmake_source=origin,backup=str(backup) if backup else None,
                existing_targets_preserved=True,component='capillary')
    print(json.dumps(result,indent=2))
    print('Configure: cmake -S . -B build -DPF_BUILD_SPHINXSYS_VALIDATION=OFF')
    print('Build: cmake --build build --target pf_soot_dem pf_soot_dem_test -j 6')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--soot-repo',type=Path,required=True)
    a=p.parse_args()
    try:install(a.soot_repo)
    except (ValueError,OSError) as exc:p.error(str(exc))
