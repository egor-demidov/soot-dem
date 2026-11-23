import math
import pathlib
import sys
import unittest
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'tools'))
from pf_interface import density_grid,fit_cap,reconstruct,history_summary

class InterfaceTests(unittest.TestCase):
    def test_kernel_integral_and_translation(self):
        xyz=np.array([[0.,0.,0.]])
        grid,_=density_grid(xyz,np.array([1.]),.4,.02)
        self.assertAlmostEqual(float(grid.sum())*.02**3,1.,places=4)
        shifted,_=density_grid(xyz+np.array([.2,-.1,.3]),np.array([1.]),.4,.02)
        self.assertAlmostEqual(float(shifted.sum())*.02**3,1.,places=4)

    def test_isosurface_volume_against_analytic_kernel_sphere(self):
        try: import skimage.measure
        except ImportError: self.skipTest('Optional scikit-image not installed')
        H=.4;level=.5
        factor=21/(2*math.pi*H**3)
        q=brentq(lambda q:factor*(1-q)**4*(1+4*q)-level,0,1)
        _,_,volume=reconstruct(np.zeros((1,3)),np.ones(1),H,.02)
        reference=4*math.pi*(q*H)**3/3
        self.assertLess(abs(volume/reference-1),.02)

    def test_volume_constrained_caps_with_shift(self):
        for angle in [40.,75.,120.]:
            theta=math.radians(angle);R=.8;shift=.012;dx=.025
            volume=math.pi*R**3*(2-3*math.cos(theta)+math.cos(theta)**3)/3
            phi=np.linspace(.01,theta,50)
            azi=np.linspace(0,2*math.pi,60,endpoint=False)
            p,a=np.meshgrid(phi,azi,indexing='ij')
            vertices=np.column_stack([(R*np.sin(p)*np.cos(a)).ravel()+.2,
                                      (R*np.sin(p)*np.sin(a)).ravel()-.3,
                                      (shift+R*(np.cos(p)-math.cos(theta))).ravel()])
            fit=fit_cap(vertices,volume,dx,np.array([.2,-.3]))
            self.assertEqual(fit['status'],'estimated')
            self.assertAlmostEqual(fit['angle_deg'],angle,places=5)
            self.assertAlmostEqual(fit['interface_shift'],shift,places=6)
            self.assertLess(fit['fit_rms'],1e-8)

    def test_late_angle_drift_excludes_invalid_geometry(self):
        frames=[dict(time=t,angle_deg=60+2*t,calibration_eligible=True) for t in range(13)]
        frames.append(dict(time=12,angle_deg=150,calibration_eligible=False))
        summary=history_summary(frames)
        self.assertEqual(summary['usable_frames'],13)
        self.assertEqual(summary['tail_usable_frames'],4)
        self.assertAlmostEqual(summary['tail_angle_slope_deg_per_time'],2,places=9)
        self.assertAlmostEqual(summary['tail_angle_range_deg'],6,places=9)
        self.assertFalse(summary['equilibrium_established'])

    def test_empty_exposed_cap(self):
        fit=fit_cap(np.zeros((15,3)),1.,.1,np.zeros(2))
        self.assertEqual(fit['status'],'insufficient exposed surface')

if __name__=='__main__': unittest.main()
