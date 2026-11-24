import math
from kinetics.phi_bose_portal_vacuum_depletion_envelope import K,trace_bound,solve

def test_vacuum_envelope_trace():
    m=.02;M=.1;mu=86e-9
    for s in [4*m*m,.002,.01,.1,1.,100.,1e4]:
        d=math.sqrt(max(0,(s-4*m*m)*(s-4*mu*mu)))/2
        for c in [-1.,0.,1.]:
            t=m*m+mu*mu-s/2+d*c
            tr=trace_bound(s,t)
            assert t<=1e-12 and tr>=-1e-12*s*s and tr<=(M*M-t)*s/2+1e-12
    assert math.isclose(K,3.978873577297384e-12)

def test_conditional_optimistic_depletion():
    r=solve(.02)
    assert r["survival_distributed"]>=r["survival_instant_bound"]
    assert .86<r["survival_distributed"]<.87
    assert r["pointwise_bose_cap_must_exceed_for_99pct_phi"]>7000
