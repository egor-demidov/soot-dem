import math,pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'tools'))
from young_laplace import bridge_reference
from analyze import sphere_fit
class ReferenceTests(unittest.TestCase):
 def test_bridge_similarity(self):
  a=bridge_reference(.25,.3,60,.6,.008)
  b=bridge_reference(.25*8,.3*2,60,1.2,.008)
  self.assertLess(a['residual'],1e-5)
  self.assertAlmostEqual(b['attractive_force']/a['attractive_force'],2,places=6)
  self.assertAlmostEqual(b['pressure']/a['pressure'],.5,places=6)
 def test_exact_cylindrical_bridge(self):
  R=.6; beta=math.pi/6; r=R*math.sin(beta); h=R*(1-math.cos(beta)); gap=.3; sigma=.008
  volume=2*math.pi*r*r*(gap/2+h)-2*math.pi*h*h*(R-h/3)
  ref=bridge_reference(volume,gap,60,R,sigma)
  self.assertAlmostEqual(ref['neck_radius'],r,places=6)
  self.assertAlmostEqual(ref['pressure'],sigma/r,places=6)
  self.assertAlmostEqual(ref['attractive_force'],math.pi*sigma*r,places=6)
 def test_exact_obtuse_spherical_bridge(self):
  # A spherical interface intersects the facing hemispheres and has zero axial force.
  radius=1.; gap=.5; sigma=1.; a=.6; center=radius+gap/2
  z=(center**2+a**2-radius**2)/(2*center)
  h=z-gap/2
  beta=math.asin(math.sqrt(a*a-z*z)/radius)
  theta=90+math.degrees(math.asin(z/a)-beta)
  volume=2*math.pi*(a*a*z-z**3/3-h*h*(radius-h/3))
  self.assertGreater(theta,90)
  ref=bridge_reference(volume,gap,theta,radius,sigma)
  self.assertAlmostEqual(ref['neck_radius'],a,places=7)
  self.assertAlmostEqual(ref['pressure'],2*sigma/a,places=7)
  self.assertAlmostEqual(ref['attractive_force'],0.,places=7)
 def test_obtuse_overhang_and_repulsion(self):
  ref=bridge_reference(36830.527,10.5,150,31.5,.08781)
  self.assertGreater(150+ref['filling_angle_deg'],180)
  self.assertLess(ref['attractive_force'],0)
  self.assertLess(ref['residual'],1e-6)
 def test_sphere_measurement(self):
  rows=[]
  for z in np.linspace(.15,.9,12):
   r=math.sqrt(1-z*z)
   for phi in np.linspace(0,2*math.pi,20,endpoint=False):rows.append(dict(x=r*math.cos(phi),y=r*math.sin(phi),z=z-.1,surface=1))
  fit=sphere_fit(rows,.01,'sessile')
  self.assertAlmostEqual(fit['radius_particle_centers'],1,places=9)
  self.assertAlmostEqual(fit['contact_angle_deg_particle_centers'],math.degrees(math.acos(.1)),places=9)
if __name__=='__main__':unittest.main()
