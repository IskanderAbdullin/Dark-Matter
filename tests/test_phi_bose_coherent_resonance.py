import numpy as np
from kinetics.phi_bose_coherent_resonance import band, floquet, independent, integrals, background, LAMBDAS, RHO0

def test_lame_and_independent_oscillator():
    for A in [.1,.3,1.,10.]:
        p,o,lo,w=band(A)
        assert np.isclose(lo,3*(1+A*A))
        assert w>0
        r=floquet(A,lo+w/2)
        alt=independent(A)
        assert r['mu_over_m']>0
        assert abs(alt['mu']/r['mu_over_m']-1)<1e-5
        assert abs(alt['det']-1)<1e-9
        i0,ic,_,_=integrals(A)
        assert abs(alt['w']-(2*ic/((1+A*A)*i0)-1))<1e-10
        assert floquet(A,lo+1.2*w)['disc_excess']<0

def test_adiabatic_normalization():
    for lam in LAMBDAS:
        assert abs(background(1,lam)['rho']/RHO0-1)<1e-10
        b=background(1e-11,lam)
        a2=1.001e-11
        b2=background(a2,lam)
        assert abs(np.log(b2['rho']/b['rho'])/np.log(1.001)+3*(1+b['w']))<.001
        assert 0<b['cs2']<1/3
