"""Conditional forward-loss bound allowing internal energy redistribution."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from scalar_pair_expanding_continuation_reference import e_total, moments
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_pair_selective_loss_bound_reference import loss_floor, moment_bounds, sha
from warm_scalar_pair_late_support_bound_reference import coefficient


def ratio(y,x,e0,P,extra):
    b2=(1-y)**(-2./3)
    lim=np.hypot(m,P/np.sqrt(b2))
    return extra*b2*(m+lim)*m*(e0-x)/(2*P*P*x)


def envelope(mu,e0,P,extra,tol=2e-10):
    """Return surviving fraction and the time of first support domination."""
    assert 0<=mu and e0>1. and P>0. and extra in (1.,2.)
    if ratio(0.,1.,e0,P,extra)>=1.:
        return 1./(1+mu),0.
    def switch(y,x):
        return ratio(y,x[0],e0,P,extra)-1.
    switch.terminal=True
    switch.direction=1
    def dz(y,x):
        return -mu*x*x*ratio(y,x[0],e0,P,extra)
    r=solve_ivp(dz,(0.,1.-1e-12),[1.],events=switch,rtol=tol,atol=tol*1e-3,max_step=.035)
    assert r.success and len(r.t_events[0])==1
    y=float(r.t_events[0][0])
    x=float(r.y_events[0][0][0])
    assert 0.<y<1. and 0.<x<=1. and abs(switch(y,[x]))<2e-8
    return 1./(1./x+mu*(1-y)),y


def independent_no_event(mu,e0,P,extra):
    end=1.-1e-9
    def rhs(y,x):
        return -mu*x*x*min(1.,ratio(y,x[0],e0,P,extra))
    z=solve_ivp(rhs,(0.,end),[1.],rtol=2.5e-11,atol=3e-14,max_step=.027)
    assert z.success and ratio(end,z.y[0,-1],e0,P,extra)>=1.
    return float(1./(1./z.y[0,-1]+mu*(1.-end)))


def run(shape,audit):
    s=json.loads(shape.read_text())["snapshots"][-1]
    a=json.loads(audit.read_text())["chained_endpoint_cases"][-1]
    assert s["steps"]==a["n_steps"] and abs(s["T_visible_GeV"]/a["temperature_GeV"]-1.)<1e-12
    T=s["T_visible_GeV"]
    q=np.asarray(s["lab_momentum_edges_GeV"])
    th=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    N,p2=moment_bounds(q,th)
    E=e_total(moments(q,th))
    P=q[-1]
    baseline=a["retained_and_all_particle_tail_Omega_h2"][1]
    inert=a["all_tail_Omega_h2"][1]
    retain=baseline-inert
    target=.120
    assert N>0. and E>m*N and 0.<inert<target<baseline
    assert abs(N*s["scale_factor"]**3/s["represented_species_count_a3"][1]-1.)<1e-11
    e0=E/(m*N)
    C=coefficient(P)
    rows=[]
    for g in (1.,.98):
        H=np.sqrt(np.pi**2/45)/MPL*(g*T)**2
        tau=N*C*P*P/(384*np.pi*m*m*H)
        for eta in (1.,2.):
            r=eta*p2[-1]/(2*P*P)
            def pack(F):
                mu=tau*(1.+F)**2
                x,y=envelope(mu,e0,P,eta)
                z,w=envelope(mu,e0,P,eta,tol=2.5e-11)
                independent=independent_no_event(mu,e0,P,eta)
                diff=max(abs(x/z-1.),abs(y-w),abs(x/independent-1.))
                assert diff<2e-8 and x>=1./(1+mu)-1e-12
                return dict(occupation_cap=F,positive_support_only_floor=inert+retain/(1.+mu),
                    energy_and_support_floor=inert+retain*x,moment_loss_only_floor=inert+retain*loss_floor(mu,r),
                    mu=mu,integrated_support_time_fraction_at_switch=y,relative_control=diff)
            ceiling=brentq(lambda F:pack(F)["energy_and_support_floor"]-target,0.,150.,xtol=2e-9)
            rows.append(dict(charge_relaxation=eta,photon_floor_ratio=g,zero_occupation_exposure=tau,
                occupation_cap_necessary_for_target=ceiling,controls=[pack(F) for F in (0.,.001,1.,2.,10.,40.)]))
    return dict(T_anchor_GeV=T,input_sha256=[sha(shape),sha(audit)],initial_active_E_per_stable_rest_mass=e0,
        initial_active_stable_n_GeV3=N,old_inert_tail_Omega=inert,active_upper_incident_p_GeV=P,rows=rows,
        scope="Conditional expansion bound: isolated positive-pressure active dark reservoir may redistribute energy or include inverse gains, but has neither external injection nor re-entry of old tails, stable number loss is only the specified capped isotropic pair channel, and incident support cools as P/b. Fixed photon floor and future Bose occupation cap assumed. ODE values are numerical comparison quadratures, not verified roundoff interval arithmetic or a full simulation.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shape",type=Path)
    p.add_argument("audit",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.shape,a.audit),indent=2,allow_nan=False))
