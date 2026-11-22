"""Exact microevent Riccati and symmetric bounded Strang splitting."""
import math
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import expit
from kinetics.vector_neutrino_exact_event_strang import (
    exact_extent,exact_sweep,pilot,
)
from kinetics.vector_neutrino_shell_radial_endpoints import adjusted_graph
from kinetics.vector_neutrino_reversible_event_graph import (
    synthetic_spectra,ledger,
)


def test_exact_local_riccati_matches_independent_ivp():
    g=adjusted_graph(.4,.2,4,scheme="shell_cell_integral")
    b,f,a=synthetic_spectra(g)
    for v,i,l,W in g["edges"][:6]:
        B,F,A=b[v],f[i],a[l]
        wv,wf,wa=g["wV"][v],g["wnu"][i],g["wnu"][l]
        def ode(_,y):
            x=y[0];bv=B-x/wv;fn=F+x/wf;fa=A+x/wa
            return [W*(bv*(1-fn)*(1-fa)-(1+bv)*fn*fa)]
        for dt in (.02,.5,4):
            calculated=exact_extent(B,F,A,wv,wf,wa,W,dt)
            ref=solve_ivp(ode,(0,dt),[0.],
                          rtol=2e-11,atol=2e-17).y[0,-1]
            assert math.isclose(calculated,ref,rel_tol=2e-8,abs_tol=2e-11)


def test_strang_preserves_quantum_invariants_and_entropy_no_thermostat():
    g=adjusted_graph(.4,.2,4,scheme="shell_cell_integral")
    initial=synthetic_spectra(g)
    previous=ledger(g,*initial)
    for scheme in ("forward","reverse","strang"):
        B,F,A=exact_sweep(g,*initial,dt=2,scheme=scheme)
        now=ledger(g,B,F,A)
        for key in ("energy","reduced_charge","lepton_asymmetry"):
            assert abs(now[key]-previous[key])<1e-12
        assert now["entropy"]>previous["entropy"]
        assert np.min(B)>=0
        assert 0<=np.min(F)<=np.max(F)<=1
        assert 0<=np.min(A)<=np.max(A)<=1


def test_strang_is_second_order_only_after_timestep_refinement():
    reference,_=pilot(h=.2,total_tau=4.,dt=.00390625)
    errors=[]
    for dt in (.0625,.03125,.015625):
        actual,r=pilot(h=.2,total_tau=4.,dt=dt)
        assert r["total_entropy_gain"]>0 and r["min_step_entropy_gain"]>0
        assert max(r["max_C_drift"],r["max_L_drift"],
                   r["max_E_drift"])<1e-12
        errors.append(max(float(np.max(abs(a-b))) for a,b in
                          zip(actual,reference)))
    assert errors[0]/errors[1]>2.5
    assert errors[1]/errors[2]>3.0
    assert errors[2]>1e-4


def test_exact_chemical_be_fd_fixed_point_and_empty_vector_inverse():
    g=adjusted_graph(.4,.2,4,scheme="shell_cell_integral")
    k,E=g["k"],g["E"]
    nu=expit(.12-k)
    antinu=expit(.03-k)
    b=1/np.expm1(E-.15)
    result=exact_sweep(g,b,nu,antinu,dt=4.,scheme="strang")
    assert max(float(np.max(abs(a-c))) for a,c in
               zip((b,nu,antinu),result))<2e-13
    bzero=np.zeros_like(E)
    after=exact_sweep(g,bzero,nu,antinu,dt=.1,scheme="strang")
    assert np.min(after[0])>0
    assert np.min(after[1])>=0 and np.max(after[1])<=1


def test_reject_invalid_pauli_occupation_and_negative_step():
    g=adjusted_graph(.4,.2,4,scheme="shell_cell_integral")
    b,f,a=synthetic_spectra(g)
    f[0]=1.01
    with pytest.raises(ValueError):
        exact_sweep(g,b,f,a,dt=.1)
    with pytest.raises(ValueError):
        exact_extent(.1,.4,.3,1.,1.,1.,1.,-1.)
