import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'research'))
from capillary_energy_reference import area_energy_forces,rigid_load,spherical_cap

class CapillaryEnergyTests(unittest.TestCase):
    def setUp(self):
        self.x=np.array([[.1,.2,.3],[1.,-.1,.2],[.2,1.1,-.2],[-.2,.3,1.]])
        self.tri=np.array([[0,1,2],[0,2,3],[0,3,1],[1,3,2]])
        self.coef=np.array([.008,.008,-.004,.008])
    def energy(self,x):return area_energy_forces(x,self.tri,self.coef)[0]
    def test_force_is_energy_gradient(self):
        _,force=area_energy_forces(self.x,self.tri,self.coef)
        for i in range(len(self.x)):
            for a in range(3):
                plus=self.x.copy();minus=self.x.copy();plus[i,a]+=1e-6;minus[i,a]-=1e-6
                self.assertAlmostEqual(force[i,a],-(self.energy(plus)-self.energy(minus))/2e-6,places=10)
    def test_translation_rotation_and_force_moment_identities(self):
        energy,force=area_energy_forces(self.x,self.tri,self.coef)
        self.assertLess(np.linalg.norm(force.sum(axis=0)),1e-16)
        self.assertLess(np.linalg.norm(np.cross(self.x,force).sum(axis=0)),1e-16)
        angle=.73;R=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
        changed,newforce=area_energy_forces(self.x@R.T+[3,1,-2],self.tri,self.coef)
        self.assertAlmostEqual(changed,energy,places=14)
        np.testing.assert_allclose(newforce,force@R.T,atol=1e-16)
    def test_rigid_translation_and_rotation_virtual_work(self):
        _,force=area_energy_forces(self.x,self.tri,self.coef)
        selected=np.array([1,2]);center=np.array([.3,.1,-.2]);F,T=rigid_load(self.x,force,selected,center)
        eps=1e-6
        for axis in np.eye(3):
            displacement=np.zeros_like(self.x);displacement[selected]=axis
            dE=(self.energy(self.x+eps*displacement)-self.energy(self.x-eps*displacement))/(2*eps)
            self.assertAlmostEqual(dE,-F@axis,places=10)
            displacement[selected]=np.cross(axis,self.x[selected]-center)
            dE=(self.energy(self.x+eps*displacement)-self.energy(self.x-eps*displacement))/(2*eps)
            self.assertAlmostEqual(dE,-T@axis,places=10)
    def test_capillary_uplift_converges_to_pressure_load(self):
        for theta in (60,90,120):
            errors=[]
            for rings in (12,24,48):
                R=.6;gamma=.008;x,top,disk,attached=spherical_cap(R,theta,rings,192)
                _,forces=area_energy_forces(x,top,gamma)
                uplift,_=rigid_load(x,forces,attached,[0,0,0])
                pressure_load=2*gamma/R*np.pi*(R*np.sin(np.deg2rad(theta)))**2
                errors.append(abs(uplift[2]-pressure_load)/pressure_load)
            self.assertLess(errors[-1],errors[0])
            self.assertLess(errors[-1],.04)
    def test_wetting_energy_enters_young_balance(self):
        # Circular cap parameterization at fixed volume. Young's angle is a
        # stationary point of gamma*A_LG-gamma*cos(thetaY)*A_SL.
        gamma=.008;V=.4;thetaY=np.deg2rad(60)
        def energy(theta):
            c=np.cos(theta);R=(3*V/(np.pi*(2-3*c+c**3)))**(1/3)
            return gamma*2*np.pi*R*R*(1-c)-gamma*np.cos(thetaY)*np.pi*R*R*np.sin(theta)**2
        eps=1e-5
        self.assertAlmostEqual((energy(thetaY+eps)-energy(thetaY-eps))/(2*eps),0,places=10)
        self.assertGreater(energy(thetaY+.1),energy(thetaY))
        self.assertGreater(energy(thetaY-.1),energy(thetaY))

if __name__=='__main__':unittest.main()
