"""Full 5D boosted-Pauli adaptive reference, conservative but not Γ_el."""
from kinetics.finite_mediator_5d_adaptive_reference import full_heat
from kinetics.independent_massive_chie_5d import cm_heat
import math


def test_5d_weak_heat_physical_sign_and_common_temperature_fixed_point():
    for mass in (0.,.5):
        hot_lepton=full_heat(m=mass,Tchi=.9,mV=.1,nr=5,nmu=5,nphi=12)
        hot_chi=full_heat(m=mass,Tchi=1.1,mV=.1,nr=5,nmu=5,nphi=12)
        equal=full_heat(m=mass,Tchi=1.,mV=.1,nr=5,nmu=5,nphi=12)
        assert hot_lepton>0 and hot_chi<0
        assert abs(equal)<1e-14*abs(hot_lepton)


def test_light_mediator_5d_angular_underresolution_negative_control():
    args=dict(mV=.03,nr=6,nmu=6,nphi=12,cut_kinetic_T=14.)
    ref=full_heat(**args,polar="adaptive")
    tr=full_heat(**args,polar="transfer",nang=24)
    cs=full_heat(**args,polar="cosine",nang=24)
    assert ref>0
    assert abs(tr/ref-1)>.1 and abs(cs/ref-1)>.1


def test_separate_5d_cosine_scheme_matches_existing_cm_normalization():
    M=.1;m=0.;Te=1.;Tchi=.9;G=.001;mV=.1;cut=14.
    pmax=math.sqrt((M+cut*Te)**2-M*M)
    kmax=math.sqrt((m+cut*Te)**2-m*m)
    reference=cm_heat(M,m,G,Te,Tchi,pmax=pmax,kmax=kmax,
                      nr=6,nang=6,nphi=12,mediator_mass=mV)
    new=full_heat(M=M,m=m,G=G,Te=Te,Tchi=Tchi,mV=mV,
                  cut_kinetic_T=cut,nr=6,nmu=6,nphi=12,polar="cosine",nang=6)
    assert math.isclose(new,reference,rel_tol=2e-6)
