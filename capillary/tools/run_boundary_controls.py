#!/usr/bin/env python3
"""Fresh fixed-sphere diagnostics; no force law or boundary modifications.

Defaults to two controls: cohesion with zero solid PF, and sfs/sff=3.
The supplied original fixed-sphere XML is retained as the reference geometry.
Removing solid PF removes attraction AND short-range repulsion.
"""
import argparse,csv,json,shutil,subprocess,sys
from pathlib import Path
import xml.etree.ElementTree as ET
from run_pf_dem import parse_xml


def make_control(base,sfs):
    root=ET.fromstring(ET.tostring(base))
    root.find('fluid').set('sfs',str(sfs))
    # There is no DEM contact in a one-fixed-sphere control.
    root.find('run').set('checkpoint_every_steps','1000')
    return root


def validate_base(path):
    config=parse_xml(path);root=ET.parse(path).getroot()
    if root.tag!='pf_dem' or root.find('particles') is not None:
        raise ValueError('Use the standalone <pf_dem> boundary_fixed_sphere.xml with an inline particle')
    if len(config['model']['bodies'])!=1 or not config['model']['bodies'][0]['fixed']:
        raise ValueError('Control requires exactly one fixed sphere')
    if 'soot_dem' in config['model']:raise ValueError('Control must not contain soot_dem')
    if root.find('run') is None:raise ValueError('Explicit run settings required')
    if config['model']['fluid'][4]<=0:raise ValueError('Control requires positive sff')
    return root,config


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('base',type=Path)
    ap.add_argument('--executable',required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--prepare-only',action='store_true')
    ap.add_argument('--include-baseline',action='store_true')
    a=ap.parse_args()
    try:
        root,config=validate_base(a.base.resolve())
        exe=shutil.which(a.executable) or (str(Path(a.executable).resolve()) if Path(a.executable).is_file() else None)
        if not exe and not a.prepare_only:raise ValueError('Executable not found')
        if a.output.exists():raise ValueError('Output directory already exists')
        cases=[('no_solid_pf',0.),('ratio3',3*config['model']['fluid'][4])]
        if a.include_baseline:cases.append(('baseline',config['model']['fluid'][5]))
        a.output.mkdir(parents=True);summary=[]
        runner=Path(__file__).with_name('run_pf_dem.py').resolve()
        audit=Path(__file__).with_name('audit_boundary_checkpoint.py').resolve()
        for name,sfs in cases:
            xml=a.output/(name+'.xml');case=make_control(root,sfs)
            ET.ElementTree(case).write(xml,encoding='utf-8',xml_declaration=True)
            parse_xml(xml) # fail before launching on malformed controls
            result_dir=a.output/name
            command=[sys.executable,str(runner),str(xml),'--executable',exe or a.executable,'--output',str(result_dir)]
            if a.prepare_only:command.append('--prepare-only')
            print('Case:',name,'sfs=',sfs,flush=True)
            result=subprocess.run(command)
            entry=dict(case=name,sfs=sfs,runner_exit=result.returncode)
            history=result_dir/'dem_history.csv'
            if history.is_file():
                with history.open() as f:rows=list(csv.DictReader(f))
                if rows:entry['last_record']=rows[-1]
            cp=result_dir/'dem_checkpoint.dat'
            if cp.is_file():
                entry['audit_exit']=subprocess.run([sys.executable,str(audit),str(cp),'--output',str(result_dir/'boundary_audit.json')]).returncode
            summary.append(entry)
            (a.output/'control_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('Summary:',a.output/'control_summary.json')
        return 0 if all(r['runner_exit']==0 and r.get('audit_exit',0)==0 for r in summary) else 1
    except (ValueError,OSError,ET.ParseError) as exc:ap.error(str(exc))
if __name__=='__main__':raise SystemExit(main())
