"""Pointwise Bose gain ceiling for the isolated collisionless loss-only continuation."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_pair_selective_loss_bound_reference import loss_floor, moment_bounds, sha
from warm_scalar_pair_late_support_bound_reference import coefficient, test


def angular_check():
    # Integral of 1/|p-k| over a unit ball; angles are integrated
    # analytically before an independent radial quadrature.
    k=np.array([0.,.02,.25,.8,1.,1.2,3.])
    exact=np.where(k<=1.,2*np.pi*(1.-k*k/3.),4*np.pi/(3*np.maximum(k,1e-50)))
    numeric=np.array([4*np.pi*quad(lambda x:x*x/max(x,r),0.,1.,
                                  points=[r] if 0.<r<1. else None,epsabs=1e-11)[0] for r in k])
    assert np.max(abs(numeric/exact-1.))<2e-13 and np.max(numeric)<=2*np.pi*(1+2e-13)
    return float(max(abs(numeric/exact-1.)))


def run(shape,audit):
    s=json.loads(shape.read_text())["snapshots"][-1]
    a=json.loads(audit.read_text())["chained_endpoint_cases"][-1]
    assert s["steps"]==a["n_steps"] and abs(s["T_visible_GeV"]/a["temperature_GeV"]-1.)<1e-12
    p=np.asarray(s["lab_momentum_edges_GeV"])
    P=p[-1]
    T=s["T_visible_GeV"]
    f=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    z=[np.min(t[0]+np.minimum(0.,t[1])) for t in f]
    assert z[0]>0
    F0=float(np.exp(-z[0])/(-np.expm1(-z[0])))
    S=float(expit(-z[1]))
    C=coefficient(P)
    emin=m*m/(np.hypot(m,P)+P)
    N,M=moment_bounds(p,f)
    inert=a["all_tail_Omega_h2"][1]
    active=a["retained_and_all_particle_tail_Omega_h2"][1]-inert
    ts=[]
    for g in (1.,.98):
        H=np.sqrt(np.pi**2/45)/MPL*(T*g)**2
        for spin in (16.,32.):
            # Three invariant measures, one flux factor, p2 angle
            # delta Jacobian and p1 Coulomb integral leave 32 E (2pi)^3.
            # The positive time integral over b^{-4} d ln b gives 1/4.
            rate=spin*C*P**6*S*S/(32*emin*(2*np.pi)**3*m*m)
            dose=rate/(4*H)
            q0=1./(1.+F0)
            assert dose<q0
            top=1./(q0-dose)-1.
            tau=N*C*P*P/(384*np.pi*m*m*H)*(1.+top)**2
            r=M[-1]/(2*P*P)
            check=quad(lambda u:np.exp(-4*u),0.,np.inf,epsabs=5e-12)[0]
            assert abs(check*4-1.)<2e-13 and F0<top<.01
            ts.append(dict(photon_floor_ratio=g,overcount_spin_and_label_multiplier=spin,
                integrated_pointwise_gain_coefficient_upper=dose,
                initial_Bose_peak_F0=F0,conditional_global_Bose_cap_after_feed=top,
                equal_charge_moment_and_support_stable_Omega_floor=inert+active*loss_floor(tau,r)))
    return dict(anchor_T_GeV=T,conditional_incoming_lab_support_GeV=P,
        unpolarized_parent_peak_occupation=S,minimum_final_lab_energy_GeV=emin,
        angular_radial_Coulomb_integral_defect=angular_check(),
        maximum_grid_matrix_element_over_spin_bound=test(P,C),inputs=[sha(shape),sha(audit)],controls=ts,
        scope="Pointwise positive scalar gain in homogeneous collisionless loss-only parents, no inverse, redistributing scattering, external source, old exterior scalar occupation or old tail re-entry. Counts in the old inert stable tail are added once. Newly produced scalar occupations include all outgoing modes in the gain supremum, even above the represented momentum mesh. The future photon floor and isolated setup are assumptions; not a bound on the fully coupled warm model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shape",type=Path)
    p.add_argument("audit",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.shape,a.audit),indent=2,allow_nan=False))
