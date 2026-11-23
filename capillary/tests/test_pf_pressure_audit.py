"""Independent offline pressure-audit checks; never run SPH trajectories."""
import math
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from pf_pressure_audit import particle_virials

META = dict(dx=.1, rho0=1., sound_speed=10., sff=.001)


def test_two_particle_half_virial_and_support():
    # Exactly two interacting particles and a third outside both supports.
    r=.2; H=.4; q=r/H; mass=.001
    rho=np.array([1.02,1.03,1.])
    eos,pf,pressure,got_mass,count=particle_virials([[0,0,0],[r,0,0],[2,0,0]],rho,META)
    p=100/7*(rho**7-1)
    derivative=-210/(math.pi*H**4)*q*(1-q)**3
    force=-mass**2/(rho[0]*rho[1])*(p[0]+p[1])*derivative
    pf_force=H*.001*(1-q)**2*(1+2*q-20.5*q*q)
    np.testing.assert_allclose(eos,[.5*r*force,.5*r*force,0],rtol=1e-13)
    np.testing.assert_allclose(pf,[.5*r*pf_force,.5*r*pf_force,0],rtol=1e-13)
    np.testing.assert_allclose(pressure,p)
    assert math.isclose(got_mass,mass) and count==1


def test_virial_translation_rotation_and_chunk_invariance():
    rng=np.random.default_rng(22)
    xyz=rng.uniform(-.3,.3,(37,3));rho=rng.uniform(1,1.04,37)
    reference=particle_virials(xyz,rho,META,chunk_size=1)
    rotation,_=np.linalg.qr(rng.normal(size=(3,3)))
    moved=particle_virials(xyz@rotation+[3,-2,1],rho,META,chunk_size=65536)
    for a,b in zip(reference[:3],moved[:3]):
        np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-16)
    assert reference[3:]==moved[3:]


def test_rejects_invalid_or_coincident_samples():
    with pytest.raises(ValueError,match='Coincident'):
        particle_virials([[0,0,0],[0,0,0]],[1,1],META)
    with pytest.raises(ValueError,match='density'):
        particle_virials([[0,0,0]],[0],META)
    with pytest.raises(ValueError,match='chunk'):
        particle_virials([[0,0,0]],[1],META,chunk_size=0)
