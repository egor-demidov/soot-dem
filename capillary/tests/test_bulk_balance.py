from pathlib import Path
import sys,unittest,math
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from pf_bulk_balance import balance
class BulkBalanceTests(unittest.TestCase):
 def test_zero_cohesion_and_material_scaling(self):
  self.assertEqual(balance(1e-9,1000,1000,0)['equilibrium_density_ratio'],1.)
  a=balance(1.75e-9,1258,1000,.002);b=balance(3.5e-9,1258,1000,.004)
  self.assertAlmostEqual(a['equilibrium_density_ratio'],b['equilibrium_density_ratio'])
  self.assertAlmostEqual(a['equilibrium_density_ratio'],1.14215708644,places=10)
  self.assertLess(abs(a['mechanical_bulk_pressure_pa']),1e-5)
  self.assertLess(balance(1.75e-9,1258,2000,.002)['equilibrium_density_ratio'],a['equilibrium_density_ratio'])
 def test_profile_integral_independently_by_quadrature(self):
  n=10000;integral=0.
  for i in range(n):
   q=(i+.5)/n;f=1-23.5*q*q-2*(2-23.5)*q**3+(3-23.5)*q**4
   integral+=q**3*f/n
  self.assertAlmostEqual(integral,(9-23.5)/168,places=8)
if __name__=='__main__':unittest.main()
