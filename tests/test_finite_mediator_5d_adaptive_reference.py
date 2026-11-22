"""Independent full 5D boosted-Pauli adaptive reference; not physical gamma_el."""
from kinetics.finite_mediator_5d_adaptive_reference import full_heat
from kinetics.independent_massive_chie_5d import cm_heat
import math


def test_5d_weak_heat_signs_and_exact_common_temperature_fixed_point():
    for mass in (0.,.5):
        cold=full_heat(m=mass,Tchi=.9,mV=.1,nr=5,nmu=5,nphi=12)
        hot=full_heat(m=mass,Tchi=1.1,mV=.1,nr=5,nmu=5,nphi=12)
        equal=full_heat(m=mass,Tchi=1.,mV=.1,nr=5,nmu=5,nphi=12)
        assert cold>0 and hot<0
        assert abs(equal)<1e-14*abs(cold)


def test_light_vector_full_5d_polar_underresolution_negative_control():
    args=dict(mV=.03,nr=6,nmu=6,nphi=12,cut_kinetic_T=14.)
    ref=full_heat(**args,polar="adaptive")
    tr=full_heat(**args,polar="transfer",nang=24)
    cs=full_heat(**args,polar="cosine",nang=24)
    assert ref>0
    assert abs(tr/ref-1)>.1 and abs(cs/ref-1)>.1


def test_existing_5d_cm_vs_independent_azimuth_reference_and_alias():
    M=.1;m=0.;Te=1.;Tchi=.9;G=.001;mV=.1;cut=14.
    pmax=math.sqrt((M+cut*Te)**2-M*M)
    kmax=math.sqrt((m+cut*Te)**2-m*m)

    def conventional(nphi):
        return cm_heat(M,m,G,Te,Tchi,pmax=pmax,kmax=kmax,
                       nr=6,nang=6,nphi=nphi,mediator_mass=mV)

    def independent(nphi):
        return full_heat(M=M,m=m,G=G,Te=Te,Tchi=Tchi,mV=mV,
                         cut_kinetic_T=cut,nr=6,nmu=6,nphi=nphi,
                         polar="cosine",nang=6)

    # DIFFERENT phi quadrature: previous 12-node Gaussian CM azimuth
    # is biased by ~1.58e-5 relative to its own 24-node value.
    # Do not weaken numerical accuracy assertions to conceal the bias.
    assert abs(conventional(12)/conventional(24)-1)>1e-5
    # At refined azimuth, the independent uniform-periodic average
    # and conventional Gaussian φ integration agree physically.
    assert math.isclose(independent(24),conventional(24),rel_tol=2e-6)
