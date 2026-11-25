"""Contact-elastic diagnostics: independent integral and conservative trial."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from quartic_consistency import sigma_over_m,lambda_from_sigma,gap,variance,shell_nodes
from quantum_elastic_moments import integrals
from elastic_conservative_trial import Elastic
from ideal_condensate_sound import moments

def test_quartic_conversion_and_gap():
    lam=lambda_from_sigma(.001)
    assert abs(sigma_over_m(lam)/.001-1)<1e-12
    q,w=shell_nodes(np.geomspace(1e-5,1,25),np.ones(24),4)
    m=gap(q,w,1e-9,lam)
    assert abs((m*m-1e-18)/(lam/2*variance(q,w,m))-1)<1e-9
    assert gap(q,w,1e9,lam)==1e9

def test_dilute_normalization_and_be_refinement():
    residual=[]
    for n in [40,80]:
        e=np.geomspace(1e-6,30,n+1)
        q=np.sqrt(e[:-1]*e[1:]);m=.3
        f=1/np.expm1(np.hypot(q,m)-m+.4)
        r=integrals(e,f,m,2,2)
        assert abs(r['dilute_normalization_rel_error'])<.002
        residual.append(abs(r['q2_net_over_absolute']))
    assert residual[1]<.35*residual[0],residual

def test_trial_exact_detailed_balance_and_entropy():
    m=Elastic(np.geomspace(1e-4,20,33),.1)
    eq=1/np.expm1(m.K+.3)
    assert m.vol@abs(m.rhs(eq))<1e-11
    f=np.exp(-m.K)*np.where(m.K<.2,50.,.1)
    r=m.diagnostics(f)
    assert r['number_residual']<1e-12 and r['energy_residual']<1e-12
    assert r['event_energy_error']<1e-12
    assert r['entropy_production_raw']>0
    step=m.positive_local_step(f,1.)
    assert step['success'] and abs(step['relative_number_change'])<1e-12
    assert abs(step['relative_energy_change'])<1e-12

def test_saturated_fluid_asymptotes():
    n,u,p=moments(1e-10)
    assert abs(u/(np.pi**2/30)-1)<1e-7
    assert abs(p/u-1/3)<1e-7
    n,u,p=moments(1e8)
    assert abs(p/u-2/3)<1e-6
    assert abs((u+p)/n-2.5*1.34148725725/2.61237534869)<1e-6
