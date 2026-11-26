"""Bound warm signed pair depletion without a scalar occupation cap."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp
from warm_scalar_pair_late_support_bound_reference import report as previous


def signed_identity():
    # Compare the unfactored and algebraically exact stable-loss integrands.
    rng=np.random.default_rng(1554)
    a,b=rng.uniform(0.,.5,(2,12000))
    c,d=np.exp(rng.uniform(-17.,17.,(2,12000)))
    pos=a*b*(1+c)*(1+d)
    neg=c*d*(1-a)*(1-b)
    expr=a*b*(1+c+d)+c*d*(a+b-1.)
    den=np.maximum(pos+neg,1e-300)
    delta=float(np.max(np.abs(pos-neg-expr)/den))
    assert delta<2e-15
    assert np.all(expr<=a*b*(1+c+d))
    return delta

def collision_control():
    """Independent gross pair integrals extract linear and quadratic terms."""
    from archived_scalar_pair_reaction_reference import gross,points
    from warm_scalar_pair_late_support_bound_reference import coefficient
    rows=[]
    for P in (.00582,.01746):
        E=np.hypot(m,P)
        low,hi=E-P,E+P
        p=np.unique(np.r_[np.linspace(0.,P,11),P,[.8*low,low,E-.001,E,E+.001,hi,hi*1.01]])
        def h(x):
            return .15*np.exp(-.5*((x-E)/.0026)**2)
        for cap in (1e-4,.19):
            def f(x):
                return cap*(x<E)
            vals=[]
            for v in (0.,1.,2.):
                _,a,_=gross(p,m,mp,f,lambda x:v*h(x),False,3,9,8,3)
                vals.append(2*a[1])
            _,a,_=gross(p,mp,m,h,f,True,3,9,8,3)
            Q=(vals[2]-2*vals[1]+vals[0])/2
            L=(vals[1]-vals[0])-Q
            n=points(p,mp,h,28)[2].sum()
            U=cap**2*coefficient(P)*P**4*(E/low)*n/(4*np.pi**3*m)
            assert vals[0]>0 and 0<L<U and 0<Q<2*a[1]
            rows.append(dict(stable_support_GeV=P,fermion_f=cap,
                 measured_linear_over_analytic_envelope=L/U,
                 quadratic_forward_over_inverse=Q/(2*a[1])))
    return rows


def exact_retained_fermion_max(shape):
    a,b=np.asarray(shape["entropy_variable_fermion"])
    return float(expit(-np.min(np.minimum(a,a+b))))


def anchor(shape,r,target=.120,photon_factor=1.):
    omega_phi,omega=r["retained_and_all_particle_tail_Omega_h2"]
    exterior=r["all_tail_Omega_h2"][1]
    assert shape["steps"]==r["n_steps"] and abs(shape["T_visible_GeV"]/r["temperature_GeV"]-1)<1e-12
    c=previous(shape,omega,exterior,target=target,photon_factor=photon_factor)
    P=c["support_upper_GeV"]
    C=c["spin_averaged_amplitude2_over_pCM2_bound_GeV_minus2"]
    n0=c["initial_inclusive_stable_number_GeV3"]
    ntotal=n0*(1+omega_phi*m/(omega*mp))
    E=np.hypot(P,m)
    # Stable equivalent energy lower gap at the anchor. This remains
    # strictly positive, even at zero scalar mass.
    Emin=m*m/(E+P)
    H=np.sqrt(np.pi**2/45)/MPL*(c["temperature_GeV"]*photon_factor)**2
    A=c["quadratic_loss_integral_without_Bose_enhancement"]
    # Fractional upper bound on the linear Bose-stimulated integrated
    # loss divided by the squared supremum of either stable occupancy.
    B=C*P**4*(E/Emin)*ntotal/(8*np.pi**3*m*H*n0)
    ref_occ=exact_retained_fermion_max(shape)
    vac=1/(1+A)
    critical=np.sqrt((vac-target/omega)/B) if vac>target/omega else 0.
    def bound(F):
        assert 0<F<=.5
        return dict(assumed_fermion_occupation_cap=F,
                    integrated_linear_loss_upper_fraction_of_initial=B*F*F,
                    conditional_final_stable_Omega_h2_floor=omega*max(0.,vac-B*F*F))
    return dict(temperature_GeV=c["temperature_GeV"],support_upper_GeV=P,
                photon_floor_comoving_temperature_fraction=photon_factor,
                stable_fixed_count_Omega_h2=omega,
                exact_anchor_retained_fermion_peak=ref_occ,
                initial_inclusive_stable_and_total_dark_density_GeV3=[n0,ntotal],
                vacuum_quadratic_coefficient=A,
                integrated_linear_loss_coefficient_per_F_squared=B,
                cap_threshold_above_which_Omega_0p120_is_not_excluded=critical,
                rows=[bound(F) for F in [ref_occ,1e-5,3e-5,1e-4,2e-4]])


def read(path):
    return json.loads(Path(path).read_text())


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--anchor",nargs=3,action="append",required=True,metavar=("SHAPE","COUNT_AUDIT","AUDIT_KEY"))
    p.add_argument("--photon-factor",type=float,default=1.)
    p.add_argument("--verify-collision",action="store_true")
    a=p.parse_args()
    results=[anchor(read(s)["snapshots"][-1],read(q)[k][-1],
                    photon_factor=a.photon_factor) for s,q,k in a.anchor]
    print(json.dumps(dict(signed_identity_relative_defect=signed_identity(),anchors=results,
          independent_collision_controls=collision_control() if a.verify_collision else [],
          scope="Prescribed late expansion floor, stable incident support inside P/b, uniform stable FD cap <=1/2, conserved total dark number and only signed stable pair conversion changing stable count. No bound on scalar occupancy or scalar momentum support. Growing stable support or additional number-changing channels can invalidate this result."),
          indent=2,allow_nan=False))
