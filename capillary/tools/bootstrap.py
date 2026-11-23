#!/usr/bin/env python3
"""Fetch only the reviewed upstream revision; never overwrite a dirty checkout."""
import pathlib, subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
PIN='319520e02f631847f85160c0e8b55648c395eafa'
p=ROOT/'deps/SPHinXsys'
def run(*args): return subprocess.check_output(args,text=True).strip()
if not p.exists():
 p.parent.mkdir(parents=True,exist_ok=True)
 subprocess.run(['git','init',str(p)],check=True)
 subprocess.run(['git','-C',str(p),'remote','add','origin','https://github.com/Xiangyu-Hu/SPHinXsys.git'],check=True)
if not (p/'.git').exists(): raise SystemExit('Dependency directory exists but is not a Git checkout')
if run('git','-C',str(p),'status','--porcelain'): raise SystemExit('Dirty dependency: refusing checkout')
subprocess.run(['git','-C',str(p),'fetch','--depth','1','origin',PIN],check=True)
subprocess.run(['git','-C',str(p),'checkout','--detach',PIN],check=True)
print('SPHinXsys pinned to',run('git','-C',str(p),'rev-parse','HEAD'))
