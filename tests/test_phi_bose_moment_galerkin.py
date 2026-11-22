"""Independent physical Bose infrared integration and entropy-dual collision tests."""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import zeta
from kinetics.phi_bose_moment_galerkin import (
    Mode,physical_moments,yukawa_moment_source,yukawa_kinematics,
)


def reference_thermal_number(m=86e-9,T=1,upper=12.1):
    pmax=math.sqrt(upper**2-m*m)
    return quad(lambda p:p*p/math.expm1(math.hypot(p,m)/T),0,pmax,
                epsabs=3e-13,epsrel=2e-12,limit=400)[0]/(2*math.pi**2)


def test_86ev_physical_bose_number_is_thermal_on_three_grids_without_fake_BEC():
    m=86e-9
    ref=reference_thermal_number(m)
    for n in (16,64,256):
        grid=np.linspace(m,12.1,n+1)
        out=physical_moments(Mode(m,grid),grid,order=36)
        assert abs(out['N']/ref-1)<1e-8
        assert abs(out['U']-sum(grid*out['M']))<1e-13
        assert np.all(out['M']>=0)
        assert math.isfinite(out['S'])
        assert out['cond_number']>0


def test_massless_bose_is_finite_and_thermal_zeta_three():
    for n in (16,64):
        grid=np.linspace(0,50,n+1)
        out=physical_moments(Mode(0,grid),grid,order=48)
        assert abs(out['N']/(zeta(3)/math.pi**2)-1)<5e-8
        assert math.isfinite(out['S']) and out['S']>0


def test_dual_jacobian_is_symmetric_negative_and_independent_energy_moments():
    for fermion,m in ((False,86e-9),(True,.1)):
        grid=np.linspace(m,8,num=25)
        eta=grid/1.1+(.02 if fermion else 0.)
        out=physical_moments(Mode(m,grid,fermion=fermion),eta,order=36)
        J=out['dM_deta']
        assert np.max(abs(J-J.T))<1e-12
        assert np.max(np.linalg.eigvalsh(J))<0
        assert abs(out['N']-out['M'].sum())<1e-15
        assert abs(out['U']-np.dot(grid,out['M']))<1e-15


def test_three_species_exact_chemical_fixed_point_and_dynamic_entropy_identity():
    m1,m2,mphi=.1,.02,86e-9
    n=64
    modes=[Mode(m,np.linspace(m,12.1,n+1),fermion=(i<2),degeneracy=1)
           for i,m in enumerate((m1,m2,mphi))]
    events=[(p,c,2e-6) for p in (0.,.03,.3,1.0) for c in (-.85,-.1,.5,.95)]
    # Thermally valid closed 1<->2: mu1 = mu2 + mu_phi, all fugacities physical.
    mus=(.005,.012,-.007)
    eta=[(mo.nodes-mu)/.8 for mo,mu in zip(modes,mus)]
    fixed=yukawa_moment_source(modes,eta,events)
    assert max(abs(v) for v in fixed['dN'])<3e-18
    assert abs(fixed['entropy_from_event'])<1e-18
    # Slightly nonthermal mother; sign can be either but entropy >=0.
    eta[0]=eta[0]+.3*np.sin(3*modes[0].nodes)
    run=yukawa_moment_source(modes,eta,events)
    assert abs(run['dark_fermion_charge'])<1e-17
    assert abs(run['parent_plus_phi_charge'])<1e-17
    assert abs(run['total_energy_residual'])<1e-16
    assert run['entropy_from_event']>=0
    assert abs(run['entropy_from_dual_gradient']-run['entropy_from_event'])<1e-16


def test_massless_yukawa_photon_like_phi_no_node_divergence_in_thermal_moments():
    m1,m2=.1,.02
    nodes1=np.linspace(m1,5,17)
    nodes2=np.linspace(m2,5,17)
    nodesphi=np.linspace(0.,5,17)
    modes=[Mode(m1,nodes1,fermion=True),Mode(m2,nodes2,fermion=True),
           Mode(0.,nodesphi)]
    eta=[m.nodes for m in modes]
    p=0.03
    e1,e2,eph=yukawa_kinematics(p,.25,m1,m2,0.)
    assert abs(e1-e2-eph)<1e-14
    out=yukawa_moment_source(modes,eta,[(p,.25,1e-6)])
    assert abs(out['rate_phi'])<1e-16


def test_dropping_first_cell_thermal_average_into_nodal_logfugacity_breaks_equilibrium():
    from kinetics.phi_bose_moment_galerkin import audit
    x=audit(bins=(16,64,256))
    assert abs(x['consistent_thermal_Rphi'])<1e-12
    assert abs(x['cell_average_substituted_Rphi'])>1e-3
    assert x['cell_average_substituted_entropy']>0
    assert x['rows'][0]['bad_rest_node_relative_factor']>1000
    assert x['rows'][0]['Galerkin_relative_error']<1e-8


def test_thermal_bose_energy_matches_independent_physical_p_integral():
    m=86e-9;T=1.;upper=12.1
    pmax=math.sqrt(upper**2-m*m)
    reference=quad(lambda p:p*p*math.hypot(p,m)/math.expm1(math.hypot(p,m)/T),
                   0,pmax,epsabs=3e-13,epsrel=2e-12,limit=400)[0]/(2*math.pi**2)
    for n in (16,64,256):
        mode=Mode(m,np.linspace(m,upper,n+1))
        out=physical_moments(mode,mode.nodes/T,order=32)
        assert abs(out['U']/reference-1)<1e-9


def test_symmetric_negative_moment_jacobian_matches_finite_difference():
    mass=.001
    mode=Mode(mass,np.linspace(mass,4.,num=13))
    eta=mode.nodes/1.1 + .08*np.sin(2*mode.nodes)
    ref=physical_moments(mode,eta,order=32)
    for j in (0,3,7,12):
        eps=1e-5
        plus=eta.copy();plus[j]+=eps
        minus=eta.copy();minus[j]-=eps
        # Eta0 initially strictly >0 thanks to m>0 and sin contribution.
        assert minus[0]>0
        d=(physical_moments(mode,plus,order=32)['M']-
           physical_moments(mode,minus,order=32)['M'])/(2*eps)
        assert np.max(abs(d-ref['dM_deta'][:,j]))<3e-9
