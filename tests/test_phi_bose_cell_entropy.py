"""Independent checks of cellwise physical Bose p-space moments and Yukawa H."""
import math
import numpy as np
import pytest
from scipy.integrate import quad
from kinetics.phi_bose_cell_entropy import (
    CellMode,cell_hats,mode_thermal_eta,mode_moments,source,
    parent_events,fixture,tangent_audit,implicit_step,ledger)


def test_cell_hats_sum_one_and_reproduce_onshell_energy():
    x=CellMode(86e-9,np.r_[86e-9,.001,.02,.4,2.2])
    for en in (86e-9,.00002,.015,.1,1.,2.2):
        k,a,b=cell_hats(x,en)
        assert 0<=a<=1 and 0<=b<=1
        assert a+b==pytest.approx(1,abs=1e-14)
        assert a*x.edges[k//2]+b*x.edges[k//2+1]==pytest.approx(en,abs=1e-13)


def test_bose_thermal_particle_number_is_physical_no_false_ir_cell():
    for n in (8,16,32):
        m=.000000086
        edges=m+(2.2-m)*np.linspace(0,1,n+1)**2
        mode=CellMode(m,edges)
        eta=mode_thermal_eta(mode)
        result=mode_moments(mode,eta,order=32)
        pmax=math.sqrt(2.2**2-m*m)
        expected=quad(lambda p:p*p/math.expm1(math.hypot(p,m)),0,pmax,epsabs=5e-12)[0]/(2*math.pi**2)
        assert result['N']==pytest.approx(expected,rel=5e-12)
        assert result['S']>0 and np.all(result['M']>=0)
        J=result['J']
        assert np.max(abs(J-J.T))<1e-17
        assert np.all(np.linalg.eigvalsh(J)<0)
        assert np.max(abs(J[0:2,2:]))==0


def test_common_T_equilibrium_nonzero_chemical_affinity_zero():
    modes,eta,events=fixture()
    for i,m in enumerate(modes):
        mu=(.03,.05,-.02)[i]
        eta[i]=mode_thermal_eta(m,T=.8,mu=mu)
    C,_,H=source(modes,eta,events)
    assert max(abs(C))<1e-25
    assert abs(H)<1e-36


def test_physical_moving_yukawa_conservation_and_micro_H():
    for n in (8,16,32):
        modes,eta,events=fixture(n=n)
        C,JC,H=source(modes,eta,events,derivative=True)
        M=[mode_moments(m,z) for m,z in zip(modes,eta)]
        s=np.cumsum([0]+[m.size for m in modes])
        dN=[sum(C[s[i]:s[i+1]]) for i in range(3)]
        dU=[np.dot(M[i]['node_energy'],C[s[i]:s[i+1]]) for i in range(3)]
        assert abs(dN[0]+dN[1])<1e-22
        assert abs(dN[0]+dN[2])<1e-22
        assert abs(sum(dU))<1e-22
        assert H>0 and sum(np.dot(z,C[s[i]:s[i+1]]) for i,z in enumerate(eta))==pytest.approx(H,rel=5e-12)
        assert np.all(np.isfinite(JC))


def test_first_cell_zero_source_stays_zero_tangent_and_implicit_dt30():
    x=tangent_audit(n=16)
    assert x['Cphi_first_cell']==[0.0,0.0]
    assert x['eta_rate_first_cell']==[-0.0,-0.0]
    assert x['phi_net_rate']<0 and x['H']>0
    modes,eta,events=fixture()
    for dt in (20.,30.,1000.):
        res=implicit_step(modes,eta,events,dt,order=24)
        assert np.allclose(res['eta'][2][:2],eta[2][:2],atol=1e-12,rtol=0)
        assert res['delta']['S']>0 and res['instant_H']>0
        assert res['delta']['S']>=res['entropy_lower_bound']-2e-10
        assert max(abs(res['delta'][x]) for x in ('N12','N1phi','U'))<2e-13
        assert min(res['eta'][2])>0
        assert res['delta']['Nphi']<0


def test_not_all_dg_grids_have_no_first_cell_source():
    x=tangent_audit(n=8)
    assert x['phi_net_rate']<0
    assert x['Cphi_first_cell'][0]<0
    assert x['Cphi_first_cell'][1]<0
    # No blanket assertion that all short dt or alternative meshes are safe.


def test_linearized_susceptibility_2x2_offdiagonal_zero_intercell():
    modes,eta,events=fixture()
    x=mode_moments(modes[2],eta[2],order=24)
    C,_,_=source(modes,eta,events)
    n=modes[2].size
    rate=np.linalg.solve(x['J'],C[-n:])
    assert np.max(abs(rate[:2]))<1e-20
    eps=.01
    yp=eta[2]+eps*rate
    nxt=mode_moments(modes[2],yp,order=24)
    assert np.max(abs(nxt['M']-x['M']-eps*C[-n:]))<1e-16
