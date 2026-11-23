from pathlib import Path
import sys,math,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from pf_wall_balance import wall_balance
class WallBalanceTests(unittest.TestCase):
 def test_bulk_balanced_wall_threshold(self):
  a=wall_balance(1.75e-9,1258,1000,.002,.006)
  b=wall_balance(1.75e-9,1258,1000,.002,.012)
  self.assertGreater(a['net_normal_N'],0);self.assertLess(b['net_normal_N'],0)
  self.assertAlmostEqual(a['critical_sfs_over_sff'],(71/16)*a['rho_ratio'],places=12)
 def test_independent_radial_integrals(self):
  # Integrate the actual polynomials and Wendland derivative at a half space.
  n=20000;ff=fs=grad=0.
  for i in range(n):
   q=(i+.5)/n
   for alpha,target in [(-23.5,'ff'),(-11.,'fs')]:
    f=1+alpha*q*q-2*(2+alpha)*q**3+(3+alpha)*q**4
    if target=='ff':ff+=q*q*f/n
    else:fs+=q*q*f/n
   grad+=210*q**3*(1-q)**3/n
  self.assertAlmostEqual(ff,-13.5/105,places=8)
  self.assertAlmostEqual(fs,-1/105,places=8)
  self.assertAlmostEqual(grad,1.5,places=8)
if __name__=='__main__':unittest.main()
