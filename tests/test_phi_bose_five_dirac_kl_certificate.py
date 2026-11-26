"""Independent five-population four-Yukawa quantum KL exact identity."""
import math
import numpy as np
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import (
    FullDiracPhysicalMoments
)
from kinetics.phi_bose_five_dirac_kl_certificate import quantum_KL_certificate


def _m():
    return FullDiracPhysicalMoments(
        P=2.,Ppair=2.,N1=3,N2=3,Nphi=8,nq=18,nrad=4,nangle=8,
        pair_nrad=4,pair_nrelative=8,pair_ncm=32,pair_naz=4
    )


def test_physically_independent_five_species_KL_plus_offroot_residual_work():
    m=_m()
    old=m.initial_etas()
    trial=[e*(1+.02*np.sin(.11*np.arange(len(e))+.5)) for e in old]
    x=quantum_KL_certificate(m,old,trial,2500.)
    assert len(x["physical_quantum_KL_by_species"])==5
    assert all(z>=-1e-11 for z in x["physical_quantum_KL_by_species"])
    assert sum(x["physical_quantum_KL_by_species"])>1e-10
    assert abs(x["root_residual_work"])>1e-12
    assert abs(x["exact_entropy_decomposition_defect"])<3e-10
    assert abs(x["reconstructed_certificate"]-
               x["entropy_certificate"])<3e-10


def test_conditional_implicit_five_species_root_entropy_gap_and_conservation():
    m=_m()
    old=m.initial_etas()
    step=m.adaptive_step(old,5000.)
    x=quantum_KL_certificate(m,old,step["etas"],step["accepted_dt"])
    assert len(x["physical_quantum_KL_by_species"])==5
    assert all(z>=-1e-11 for z in x["physical_quantum_KL_by_species"])
    assert x["microscopic_H_new"]>0
    assert x["entropy_gain"]>0
    assert x["entropy_certificate"]>=-1e-10
    assert abs(x["exact_entropy_decomposition_defect"])<1e-10
    assert math.isclose(x["entropy_certificate"],
                        step["entropy_certificate"],abs_tol=1e-10,rel_tol=2e-6)
    assert abs(x["U1_defect"])<3e-11
    assert abs(x["energy_defect"])<3e-11


def test_full_U1_thermal_equilibrium_zero_five_species_KL():
    m=_m()
    old=m.initial_etas(T=(1.,)*5,mu=(.05,-.05,.05,-.05,0.))
    out=m.implicit_step(old,1e4)
    x=quantum_KL_certificate(m,old,out["etas"],1e4)
    assert all(abs(z)<1e-17 for z in x["physical_quantum_KL_by_species"])
    assert abs(x["entropy_gain"])<1e-20
    assert abs(x["microscopic_H_new"])<1e-26
    assert abs(x["root_residual_work"])<1e-17
