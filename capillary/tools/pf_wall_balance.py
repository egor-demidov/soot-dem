"""Static planar-wall diagnostic for the published dummy-particle PF model.

Assumes a sharp plane, homogeneous exterior liquid, fixed rho0 ghost samples,
zero velocity, continuum quadrature, and no gas. This is NOT an exact discrete
sphere penetration criterion or a contact-angle calibration.
"""
import argparse,math,json
from pf_bulk_balance import balance

def wall_balance(dx,rho0,sound,sff,sfs,density_ratio=None):
 bulk=balance(dx,rho0,sound,sff)
 if not math.isfinite(sfs) or sfs<0:raise ValueError('Invalid sfs')
 r=bulk['equilibrium_density_ratio'] if density_ratio is None else density_ratio
 if not math.isfinite(r) or r<=0:raise ValueError('Invalid density ratio')
 P=rho0*sound*sound/7*(r**7-1)
 # integral_0^1 q^2 f(q) dq = (10+alpha)/105.
 # grad kernel over a half space gives 3/(2H). Forces are on ONE liquid sample.
 pressure=-3/8*P*dx*dx/r
 liquid_pair=math.pi*4**4*13.5/105*sff*dx*r
 solid_pair=-math.pi*4**4/105*sfs*dx
 coeff=math.pi*4**4/105*dx
 critical_sfs=(pressure+liquid_pair)/coeff
 return dict(assumptions='continuum planar uniform liquid and fixed rho0 dummy solid; zero velocity',rho_ratio=r,
             pressure_normal_N=pressure,liquid_pair_normal_N=liquid_pair,solid_pair_normal_N=solid_pair,
             net_normal_N=pressure+liquid_pair+solid_pair,
             critical_sfs=critical_sfs,critical_sfs_over_sff=(critical_sfs/sff if sff else None),
             supplied_sfs_over_sff=(sfs/sff if sff else None),normal_sign='positive outward from solid')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name,default in [('dx',1.75e-9),('rho0',1258),('sound',1000),('sff',.002),('sfs',.012)]:p.add_argument('--'+name,type=float,default=default)
 p.add_argument('--density-ratio',type=float)
 a=p.parse_args();print(json.dumps(wall_balance(a.dx,a.rho0,a.sound,a.sff,a.sfs,a.density_ratio),indent=2))
