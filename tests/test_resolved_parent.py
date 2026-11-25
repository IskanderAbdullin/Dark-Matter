"""Momentum-resolved detailed balance, conserved events and spin sensitivity."""
import sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'kinetics'))
from resolved_parent_bath import ResolvedParent
from yukawa_rate_bath import YukawaRateBath
from laine_portal_rates import PortalTable,make_models

@pytest.mark.parametrize('g',[2,4])
def test_resolved_operator(g):
    model=ResolvedParent(n=64,nparent=64,g_fermion=g,kappa=.00527818,
        scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8)
    s=model.validate()
    assert s['balance']<1e-10 and s['energy']<1e-12 and s['number']<1e-12
    # Discrete H-theorem for exact positive forward/backward fluxes.
    model.geometry(.2)
    rng=np.random.default_rng(192)
    for i in range(3):
        p=rng.uniform(.02,10,model.n);f=rng.uniform(.01,.8,model.n+model.n1)
        y=np.r_[p,f]
        rhs=model.dark(y)
        grad=np.r_[np.log1p(1/p),np.log((1-f)/f)]
        ent=(model.vols*grad)@rhs
        assert ent>=0

def test_two_spin_number_closure():
    m=YukawaRateBath(n=48,g_fermion=2,kappa=1e-5,kappa_sm1=47.,
                    scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8,kappa2=.03)
    r=m.run(.1,.3,2e-6,1e-15,True,86000.,'empty')
    assert r['g_fermion']==2 and r['frozen_cell_spectrum']['g_chi2']==2
    assert abs(r['final']['continuity_residual'])<1e-5
    assert abs(r['final']['number_source_residual'])<1e-7

@pytest.mark.skipif(not (ROOT/'.runtime/mesh_imPiR.dat').exists(),reason='External reference mesh not redistributed')
def test_published_rate_convention_and_jacobian():
    t=PortalTable(ROOT/'.runtime/mesh_imPiR.dat')
    q=np.array([.1,.2,1.,3.,10.])
    assert np.max(abs(t.factor(1e4,40.,q)-1))<1e-12
    f=t.factor(860000.,.1,q)
    assert np.all(f>0) and np.ptp(f)>1
    C,R=make_models()
    kw=dict(table=t,n=64,scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8,
            kappa=.00527818,kappa_sm1=47.,kappa2=.03,g_fermion=2)
    assert C(**kw).test_full()['jacobian']<1e-5
    assert R(nparent=64,**kw).validate()['jacobian']<1e-5


def test_handoff_gate_rejects_early_closure():
    import json
    from hybrid_parent_closure import handoff_passes
    data=json.loads((ROOT/'kinetics/results/laine-and-resolved-parent-2026-11-25.json').read_text())
    assert not handoff_passes(data['hybrid64_early.json'][0]['handoff'])
    assert handoff_passes(data['hybrid120_reg.json'][0]['handoff'])


def test_identical_majorana_pair_multiplicity():
    # Two independent identical pairs share the same scalar distribution.
    # Dark production is doubled for the scalar but each fermion occupation
    # obeys the same equation. The background must use the total density.
    kw=dict(n=64,nparent=64,kappa=.006,scalar_ratio=1e-12,
            fermion_ratio=.9,qmin=1e-8,kappa_sm1=47.,kappa2=.03)
    m2=ResolvedParent(g_fermion=2,**kw);m4=ResolvedParent(g_fermion=4,**kw)
    rng=np.random.default_rng(61);n=m2.n
    y=np.r_[rng.uniform(.01,3.,n),rng.uniform(.001,.4,n+m2.n1)]
    for x in [.1,.3,2.,10.]:
        m2.geometry(x);m4.geometry(x)
        c2=m2.dark(y);c4=m4.dark(y)
        assert np.allclose(c4[:n],2*c2[:n],rtol=1e-13,atol=1e-15)
        assert np.allclose(c4[n:],c2[n:],rtol=1e-13,atol=1e-15)
        assert np.allclose(m4.energy_weights[n:],2*m2.energy_weights[n:])
        assert np.allclose(m4.energy_weights[:n],m2.energy_weights[:n])
        assert np.allclose(m4.rates(x,1.,.99)[0],m2.rates(x,1.,.99)[0])
