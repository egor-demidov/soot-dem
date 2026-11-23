"""Static continuum estimate of bulk PF cohesion balanced by the Tait EOS.
Not a droplet simulation, contact-angle fit, or exact finite-lattice prediction.
"""
import argparse,math,json

def balance(dx,rho0,sound,sff):
 if not all(math.isfinite(x) and x>0 for x in (dx,rho0,sound)) or not math.isfinite(sff) or sff<0:raise ValueError('Invalid material input')
 # Integral_0^1 q^3 f_ff(q) dq = (9+alpha_ff)/168.
 C=-(2*math.pi/3)*4**5*((9-23.5)/168)
 A=C*sff/dx
 def eos(r):return sound*sound*rho0/7*(r**7-1)
 lo,hi=1.,2.
 while eos(hi)<A*hi*hi:hi*=2
 for _ in range(100):
  mid=(lo+hi)/2
  if eos(mid)<A*mid*mid:lo=mid
  else:hi=mid
 r=(lo+hi)/2
 return dict(cohesive_coefficient=C,reference_cohesive_stress_pa=-A,equilibrium_density_ratio=r,equilibrium_eos_pressure_pa=eos(r),equilibrium_cohesive_pressure_pa=-A*r*r,mechanical_bulk_pressure_pa=eos(r)-A*r*r)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--dx',type=float,default=1.75e-9);p.add_argument('--rho0',type=float,default=1258)
 p.add_argument('--sound',type=float,default=1000);p.add_argument('--sff',type=float,default=.002)
 a=p.parse_args();print(json.dumps(balance(a.dx,a.rho0,a.sound,a.sff),indent=2))
