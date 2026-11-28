"""Enclose the conditional homogeneous expansion difference without subtracting H."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq
from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp
from warm_light_occupation_stage_audit_reference import read


def run(eos,psd,metric):
    h=.6736
    omegaD=.12/h**2
    bg=Background(eos,omega_m=(.12+.02237)/h**2)
    R=read(psd)
    beta=read(metric)["held_and_actual_direct_scalar_rest_fractions_with_same_archived_stable_count"][1]
    e=np.array([x if x is not None else 1300. for x in R["conserved_momentum_edges_q_over_T0"]])
    # Envelope for the retained Q<=max(60T,V+30T)+75T finite
    # lepton quadrature from 12–50 MeV; pion measured support is smaller.
    assert (10+105*.05)/.012<e[-1]
    Y=np.asarray(R["controls"][-1]["shell_yields"]).sum(axis=0)
    Y/=Y.sum()
    q=e*bg.T0/mp
    def kinetic(a):
        a=np.asarray(a)
        x=q[None,:]/a[:,None]
        cur=x*x/(np.hypot(1.,x)+1.)
        today=q*q/(np.hypot(1.,q)+1.)
        # h_ncdm(1) normalization is bounded independently of kinetic
        # excess at a. The constituent is replacing baseline CDM, not
        # being added on top of it.
        lo=np.maximum(0.,cur[:,:-1]@Y-Y@today[1:])/(1+Y@today[1:])
        hi=(cur[:,1:]@Y-Y@today[:-1])/(1+Y@today[:-1])
        at=a==1.
        lo[at]=0.;hi[at]=0.
        return np.stack((lo,hi),axis=1)
    def stable_delta(a):
        lo,hi=kinetic(a).T
        d=omegaD*beta*a
        coeff=bg.weighted_radiation(a)+bg.om*a+bg.ol*a**4
        L=d*lo/coeff;H=d*hi/coeff
        return np.stack((L/(np.sqrt(1+L)+1),H/(np.sqrt(1+H)+1)),axis=1)
    a0=bg.scale(.05)
    sample=np.r_[a0,bg.scale(.012),bg.scale(.0073940494651976),1e-7,1e-6,1e-5,bg.aeq,1/21,1.]
    D=stable_delta(sample)
    # Equality compares baseline present-normalized matter with
    # radiation plus the scalar's energy excess over its CDM replacement.
    shifts=[]
    for side in (0,1):
        def eq(d):
            a=bg.aeq*(1+d)
            return d-(omegaD*beta/bg.om)*(1+d)*kinetic(np.array([a]))[0,side]
        shifts.append(float(brentq(eq,0.,1e-5,xtol=1e-18)))
    assert 0.<shifts[0]<shifts[1]
    today=q*q/(np.hypot(1.,q)+1.)
    neff=(omegaD*beta/bg.ogamma)/(7/8*(4/11)**(4/3))*np.array([Y@q[:-1]/(1+Y@today[1:]),Y@q[1:]/(1+Y@today[:-1])])
    time=[]
    for order in (6,12):
        bounds=np.unique(np.r_[np.log([a0,bg.aeq]),np.clip(bg.v,*np.log([a0,bg.aeq])),
                               np.linspace(*np.log([a0,bg.aeq]),40)])
        x,w=roots_legendre(order)
        half=np.diff(bounds)[:,None]/2
        t=(bounds[:-1,None]+half*(1+x)).ravel()
        V=np.exp(t)
        w=(half*w).ravel()*V
        denom=np.sqrt(bg.weighted_radiation(V)+bg.om*V+bg.ol*V**4)
        ref=float(bg.D0*(w@(1/denom)))
        H=stable_delta(V)
        # Exact stable differences 1/(1+d)-1 = -d/(1+d).
        alt=-bg.D0*((w/denom)@(H/(1.+H)))
        time.append(dict(quadrature_order=order,baseline_reference_conformal_duration_Mpc=ref,
                         conditional_conformal_time_difference_bounds_Mpc=[float(alt[1]),float(alt[0])]))
    return dict(scope="Positive-shell homogeneous energy envelope for the conditional direct scalar fraction "
                      "replacing CDM in the same fixed background, with the scalar held collisionless. "
                      "Counterfactual treat spectrum as already deposited at earliest reference a; no actual early "
                      "production history, fixed-perturbation response or inference of dark abundance.",
                physical_rest_fraction=beta,conservative_overflow_q_over_today_T=1300.,
                scale_factors=sample.tolist(),conditional_relative_Hubble_difference_bounds=D.tolist(),
                reference_equality_relative_scale_shift_bounds=shifts,relativistic_limit_Delta_N_eff_bounds=neff.tolist(),
                conformal_time_controls=time)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","psd","metric"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.psd,a.metric),indent=2,allow_nan=False))
