"""Independent finite-mass five-dimensional CM normalization/Pauli smoke tests."""
import math
from kinetics.independent_massive_chie_5d import cm_heat, transfer_heat

def test_massless_5d_matches_prior_independent_reference():
    # Old independently measured continuum Gauss n28 benchmark with 22 GeV cutoff.
    Q = cm_heat(.1,0.,.001,1.,.9,pmax=22.,kmax=22.,nr=18,nang=18,nphi=18)
    assert math.isclose(Q,1.9851492233308662e-8,rel_tol=5e-5)

def test_massive_5d_against_exact_invariant_q_not_itself():
    m=.00051099895; M=.1*m;T=m
    pmax=math.sqrt((M+30*T)**2-M*M)
    kmax=math.sqrt((m+30*T)**2-m*m)
    q5=cm_heat(M,m,1.,T,.9*T,pmax=pmax,kmax=kmax,nr=18,nang=18,nphi=18)
    qe=transfer_heat(M,m,1.,T,.9*T,n=36,nq=16)
    assert q5>0 and qe>0 and abs(qe/q5-1)<.003
    qhot=cm_heat(M,m,1.,T,1.1*T,pmax=pmax,kmax=kmax,
                 nr=18,nang=18,nphi=18)
    assert qhot<0

def test_massive_equal_temperature_no_fake_heat():
    m=.00051099895;T=m;M=.1*m
    pmax=math.sqrt((M+25*T)**2-M*M)
    kmax=math.sqrt((m+25*T)**2-m*m)
    # Roundoff relative to typical off-equilibrium ~4e-32 GeV^5.
    assert abs(cm_heat(M,m,1.,T,T,pmax=pmax,kmax=kmax,
                       nr=12,nang=12,nphi=12))<1e-44
    assert abs(transfer_heat(M,m,1.,T,T,n=16,nq=12))<1e-44
