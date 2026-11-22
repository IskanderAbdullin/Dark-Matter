"""Negative 5D cutoff/radial convergence controls for coupled χ1χ2φℓ."""
import math
import numpy as np
from scipy.special import expit
from kinetics.four_species_yukawa_shared_lepton import configuration,initial,rhs,ledger,guarded_step
from kinetics.chi1_chi2_shared_lepton_pauli import interp_fermi
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def test_segmented_elastic_heat_matches_independent_5d_same_physical_f():
    c=configuration(bins=24,elastic_radial_scheme='segmented_p',elastic_segments=6,
        elastic_nr=4,elastic_nmu=8,elastic_nc=16,elastic_nphi=16)
    f=initial(c,T1=1.1,T2=1.,Tlep=1.)
    d,o=rhs(c,f)
    ff=interp_fermi(c['elastic']['chi_grid'][0],f[0])
    ll=interp_fermi(c['elastic']['lg'],f[2])
    independent=evaluate((ff,ll),M=.1,m=0.,g=.01,mv=.1,
       nr=4,nmu=8,nc=16,nphi=16,cutoff=6,
       radial_scheme='segmented_p',segments=6,angular_scheme='log_transfer')
    assert math.isclose(o['elastic_Qchi'][0],independent['Qchi'],rel_tol=3e-11)
    assert math.isclose(o['elastic_entropy'],independent['entropy'],rel_tol=3e-11)
    assert o['entropy_identity_abs_error']<3e-24
    assert max(abs(z) for z in o['charge_residuals'])<3e-24
    assert abs(o['energy_residual'])<3e-24
    step=guarded_step(c,f,1e6)
    assert step is not None
    assert ledger(c,step)['entropy']>ledger(c,f)['entropy']


def test_three_radial_orders_stable_locally_but_uv_cutoff_omission_large():
    rates=[]
    for n in (3,4,5):
        c=configuration(bins=24,elastic_radial_scheme='segmented_p',elastic_segments=6,
           elastic_nr=n,elastic_nmu=8,elastic_nc=16,elastic_nphi=16)
        z=rhs(c,initial(c,T1=1.1))[1]
        rates.append(z['elastic_Qchi'][0])
        assert z['elastic_entropy']>0 and z['yukawa_entropy']>0
    assert max(rates)<0 and max(abs(x/rates[-1]-1) for x in rates)<.0002
    # Original cut=6 excludes relevant high-p scattering: a distinct failure.
    extended=configuration(bins=24,elastic_cut=10.,elastic_radial_scheme='segmented_p',
         elastic_segments=10,elastic_nr=4,elastic_nmu=8,elastic_nc=16,elastic_nphi=16)
    r10=rhs(extended,initial(extended,T1=1.1))[1]['elastic_Qchi'][0]
    assert r10<rates[-1] and abs(r10/rates[-1]-1)>.20


def test_nonthermal_radial_alias_cannot_be_hidden_by_H_theorem():
    def rate(scheme,n,segments):
        c=configuration(bins=64,elastic_radial_scheme=scheme,elastic_segments=segments,
            elastic_nr=n,elastic_nmu=8,elastic_nc=16,elastic_nphi=16)
        state=list(initial(c,T1=1.1))
        state[0]=expit(-c['yukawa'].e1+.3*np.sin(4*c['yukawa'].e1))
        state[2]=expit(-c['elastic']['lg']+.3*np.cos(4*c['elastic']['lg']))
        r=rhs(c,tuple(state))[1]
        assert r['entropy_event']>0 and r['entropy_identity_abs_error']<1e-23
        return r['elastic_Qchi'][0]
    coarse=rate('uniform_p',4,1)
    refined=rate('segmented_p',4,6)
    assert coarse>2*refined>0
