"""Read-only momentum accounting. Treating CSF as external is NOT conservation."""
import csv
import numpy as np

def read(path):
    with path.open() as f: return list(csv.DictReader(f))

def vector(row, prefix):
    return np.array([float(row[prefix+a]) for a in 'xyz'])

def force_audit(folder):
    h=read(folder/'history.csv')
    if 'capillary_impulse_x' not in h[0]: return None
    p0=vector(h[0],'p')+vector(h[0],'support_impulse_')
    c0=vector(h[0],'capillary_impulse_')
    raw=[];external=[]
    for r in h:
        balance=vector(r,'p')+vector(r,'support_impulse_')-p0
        cap=vector(r,'capillary_impulse_')-c0
        raw.append(float(np.linalg.norm(balance)))
        external.append(float(np.linalg.norm(balance-cap)))
    result={'raw_max_supported_momentum_defect':max(raw),
            'max_residual_treating_CSF_as_external':max(external),
            'final_supported_momentum_defect':(vector(h[-1],'p')+vector(h[-1],'support_impulse_')-p0).tolist(),
            'final_CSF_impulse':(vector(h[-1],'capillary_impulse_')-c0).tolist(),
            'warning':'Subtracting CSF is an accounting diagnostic only. No wetting reaction or conservative force transfer is established.'}
    if 'viscous_dt_limit' in h[0]:
        ratios=[float(r['last_dt'])/float(r['viscous_dt_limit']) for r in h if float(r['viscous_dt_limit'])>0]
        result['max_recorded_dt_over_viscous_limit']=max(ratios,default=0.)
        result['viscous_limit_note']='Diagnostic at recorded times; bound applies to isolated frozen-coefficient viscosity, not all SPH terms.'
    path=folder/'momentum_stages.csv'
    if path.exists():
        groups={}
        for r in read(path): groups.setdefault(r['step'],[]).append(r)
        maxima={};local=[]
        expected=['before_pressure','after_liquid_pressure','after_gas_pressure','after_liquid_density','after_gas_density']
        for group in groups.values():
            if [r['stage'] for r in group]!=expected:
                raise ValueError('Incomplete or reordered momentum audit stages')
            def momentum(r): return vector(r,'liquid_p')+vector(r,'gas_p')
            for left,right in zip(group,group[1:]):
                kick=momentum(right)-momentum(left)
                maxima[right['stage']]=max(maxima.get(right['stage'],0),float(np.linalg.norm(kick)))
            dt=float(group[0]['dt'])
            local.append(float(np.linalg.norm(momentum(group[-1])-momentum(group[0])+
                       dt*vector(group[0],'wall_f')-dt*vector(group[0],'capillary_f'))))
        result['sampled_stage_count']=len(groups)
        result['max_stage_momentum_kicks']=maxima
        result['max_sampled_substep_residual_treating_CSF_as_external']=max(local,default=0.)
        result['stage_note']='Wall force is sampled before pressure initialization; time levels and unequal phase operators can leave a residual.'
    path=folder/'gas_outliers.csv'
    if path.exists():
        records=read(path)
        if records:
            worst=max(records,key=lambda r:float(r['speed']))
            result['fastest_recorded_gas_particle']={k:float(v) for k,v in worst.items()}
    return result
