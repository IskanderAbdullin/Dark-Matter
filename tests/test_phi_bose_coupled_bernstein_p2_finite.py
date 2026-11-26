"""Focused finite-time physical-capacity coupled χ1/χ2/φ Bernstein-P2 checks."""
import math
import numpy as np
from scipy.sparse import hstack
from kinetics.phi_bose_coupled_bernstein_p2_finite import (
    Coupled,integrate_to,conservative_adaptive_step)
from kinetics.phi_bose_physical_cutcell_weak import analytic_weak


def test_physical_modes_event_invariant_nullspace_three_grids_cutoffs():
    for P in (2.,8.,16.):
        reference=analytic_weak(P)
        for nr,na,nphi in ((4,8,8),(6,12,12),(8,16,16)):
            model=Coupled(P=P,Nphi=nphi,nrad=nr,nangle=na)
            src=model.source(model.initial_etas())
            assert np.min(model.W)>0
            assert src['H']>0
            for k in ('Rphi','Qphi','Q1','Q2'):
                assert abs(src[k]/reference[k]-1)<1e-6,(P,nphi,k)
            assert abs(src['Qphi']+src['Q1']+src['Q2'])<2e-24
            for s in model.species:
                cap_from_quadrature=np.sum(s.weights,axis=1)
                p0=np.sqrt((s.edges[:-1]-s.mass)*(s.edges[:-1]+s.mass))
                p1=np.sqrt((s.edges[1:]-s.mass)*(s.edges[1:]+s.mass))
                cap_exact=s.g*(p1**3-p0**3)/(6*math.pi**2)
                assert np.max(abs(cap_from_quadrature/cap_exact-1))<2e-13
            S=hstack([-model.A[0],model.A[1],model.A[2]]).toarray()
            sv=np.linalg.svd(S,compute_uv=False)
            rank=int(np.sum(sv>1e-11*sv[0]))
            zeros=int(np.sum(np.linalg.norm(S,axis=0)<1e-13))
            assert S.shape[1]-rank ==3+zeros,(P,nphi,S.shape[1]-rank,zeros)
            n1,n2,np0=[s.nvars for s in model.species]
            charge12=np.r_[np.ones(n1),np.ones(n2),np.zeros(np0)]
            charge1phi=np.r_[np.ones(n1),np.zeros(n2),np.ones(np0)]
            energy=np.concatenate([s.nodal_energies for s in model.species])
            assert np.max(abs(S@charge12))<1e-13
            assert np.max(abs(S@charge1phi))<1e-13
            assert np.max(abs(S@energy))<1e-12


def test_physical_finite_time_closed_three_species_charge_energy_entropy():
    for P in (2.,8.,16.):
        model=Coupled(P=P,Nphi=16,nrad=8,nangle=16)
        for variant in ('thermal_reverse','nonthermal'):
            ini=model.initial_etas()
            if variant=='nonthermal':
                ini[0]+=0.05*np.sin(ini[0]*4)
                ini[1]-=0.03*np.cos(ini[1]*2)
                ini[2]+=0.01*(1-np.cos(ini[2]))
            out=conservative_adaptive_step(model,ini,1e5)
            assert out['accepted_dt']==1e5 and len(out['rejections'])==0
            assert out['entropy_increase']>0
            assert out['source_old']['H']>0 and out['source_new']['H']>0
            assert out['min_bose_eta']>0
            assert max(abs(x) for x in out['number_residuals'])<3e-14
            assert abs(out['energy_defect'])<3e-14
            for spec,eta in zip(model.species,out['etas']):
                e=spec.evaluate(eta,with_entropy=True)
                assert np.min(e['M'])>=0
                assert math.isfinite(e['S'])
            src=model.source(ini)
            entropy_dual=sum(np.dot(eta,c) for eta,c in zip(ini,src['C']))
            assert math.isclose(entropy_dual,src['H'],rel_tol=3e-12)


def test_exact_nonzero_mu_thermal_equilibrium_on_every_event():
    model=Coupled(P=2.,Nphi=16,nrad=6,nangle=12)
    etas=model.initial_etas(T=(1.,1.,1.),mu=(-.15,.05,-.2))
    source=model.source(etas)
    assert max(abs(source[k]) for k in ('Rphi','Qphi','Q1','Q2'))<1e-24
    assert source['H']<1e-30
    assert max(np.max(abs(c)) for c in source['C'])<1e-24
    out=conservative_adaptive_step(model,etas,1e5)
    assert out['accepted_dt']==1e5
    assert max(np.max(abs(a-b)) for a,b in zip(out['etas'],etas))<1e-14


def test_adaptive_rejection_is_not_a_mathematical_strong_no_go():
    model=Coupled(P=2.,Nphi=16,nrad=8,nangle=12)
    ini=model.initial_etas()
    try:model.step(ini,3e5)
    except (ValueError,ArithmeticError):pass
    else:raise AssertionError('Intentionally oversized explicit moment step should not invert under this prototype')
    out=integrate_to(model,ini,3e5,3e5)
    assert out['accepted_steps']==2
    assert out['rejected_attempts']==1
    assert all(v['entropy_increase']>0 for v in out['history'])
    assert all(v['min_bose_eta']>0 for v in out['history'])
    assert math.isclose(out['elapsed'],3e5,rel_tol=0,abs_tol=1e-9)


def test_first_order_dt_convergence_not_false_time_independence():
    model=Coupled(P=2.,Nphi=16,nrad=8,nangle=12)
    ini=model.initial_etas()
    initial=model.integrate_spectrum(ini)
    nchanges=[]
    for maxdt in (1.875e4,9375.,4687.5):
        out=integrate_to(model,ini,3e5,maxdt)
        e=model.integrate_spectrum(out['etas'])
        nchanges.append(e['N'][2]-initial['N'][2])
    errcoarse=abs(nchanges[0]-nchanges[1])
    errfine=abs(nchanges[1]-nchanges[2])
    assert 1.9<errcoarse/errfine<2.2,(nchanges,errcoarse/errfine)
    assert errcoarse>1e-9
