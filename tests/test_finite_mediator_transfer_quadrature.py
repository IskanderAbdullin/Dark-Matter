"""Independent transfer-variable vs direct cosine quadrature of boosted Pauli heat.

This test does not claim continuum convergence or a physical cosmological rate.
"""
import math
import pytest
from kinetics.independent_massive_chie_5d import cm_heat


@pytest.mark.parametrize("mass,lepton_mass", [(1.,0.),(3.,0.),(3.,.5)])
def test_independent_transfer_angle_quadrature(mass, lepton_mass):
    M,T,G=.1,1.,.001
    pmax=math.sqrt((M+14*T)**2-M*M)
    kmax=math.sqrt((lepton_mass+14*T)**2-lepton_mass**2)
    args=dict(M=M,m=lepton_mass,G=G,Te=T,Tchi=.9*T,
              pmax=pmax,kmax=kmax,nr=12,nphi=12,mediator_mass=mass)
    a=cm_heat(**args,nang=24,angular_scheme="cosine")
    b=cm_heat(**args,nang=24,angular_scheme="transfer")
    assert math.isfinite(a) and math.isfinite(b)
    assert a>0 and b>0
    # Coarse 5D radial/azimuth grid is shared; only angular rule changes.
    assert b==pytest.approx(a,rel=.06)


def test_equal_temperature_transfer_and_input_validation():
    args=dict(M=.1,m=0.,G=.001,Te=1.,Tchi=1.,pmax=10.,kmax=10.,
              nr=10,nang=10,nphi=10,mediator_mass=.3)
    assert abs(cm_heat(**args,angular_scheme="transfer"))<1e-17
    with pytest.raises(ValueError):
        cm_heat(**args,angular_scheme="unknown")
    with pytest.raises(ValueError):
        cm_heat(**{**args,"mediator_mass":None},angular_scheme="transfer")
