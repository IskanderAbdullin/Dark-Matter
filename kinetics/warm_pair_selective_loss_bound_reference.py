"""Late shape-preserving redshift/loss-only selective stable-depletion bound."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Cells
from scalar_pair_expanding_continuation_reference import moments
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_scalar_pair_late_support_bound_reference import coefficient, test


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def loss_floor(mu,r):
    assert 0.<r<1. and mu>=0.
    transition=-np.log(r)/r
    return float(np.exp(-r*mu) if mu<=transition else 1./(1./r+mu-transition))


def mu_required(y,r):
    assert 0.<y<1.
    if y>=r:
        return float(-np.log(y)/r)
    return float(-np.log(r)/r+1./y-1./r)


def moment_bounds(q,theta):
    U=moments(q,theta)[1]
    K=q*q/(np.hypot(q,m)+m)
    N=float(U[0].sum())
    low,high=np.array([(2*m+K[:-1])@U[1],(2*m+K[1:])@U[1]])/N
    # Direct positive reconstructed f(q) cross-check, not needed for
    # the inequality. A narrow positive shell is integrated in K.
    c=Cells(q,m,4)
    e=m+K[:-1,None]+np.diff(K)[:,None]*c.u
    w=c.w*c.f(theta[1],e,np.arange(len(K)-1)[:,None])*c.vol[:,None]
    kk=K[:-1,None]+np.diff(K)[:,None]*c.u
    check=float(np.sum(w*kk*(2*m+kk))/w.sum())
    assert low<check<high
    return N,[float(low),check,float(high)]


def angle_identity(P):
    r,w=roots_legendre(16)
    x=np.geomspace(.002*P,P,25)
    A,B=np.meshgrid(x,x)
    ka=A*A/(np.hypot(m,A)+m)
    kb=B*B/(np.hypot(m,B)+m)
    em=m*(ka+kb)+ka*kb
    cm2=(em[...,None]-A[...,None]*B[...,None]*r)/2.
    assert np.all(cm2>=0.)
    bound=(A*A+B*B)/4.
    mean=cm2@(w/2)
    assert np.max(mean/bound)<1.+2e-12
    return float(np.max(mean/bound))


def run(shape,audit):
    s=json.loads(shape.read_text())["snapshots"][-1]
    d=json.loads(audit.read_text())["chained_endpoint_cases"][-1]
    assert s["steps"]==d["n_steps"] and abs(s["T_visible_GeV"]/d["temperature_GeV"]-1.)<1e-12
    total=d["retained_and_all_particle_tail_Omega_h2"][1]
    ext=d["all_tail_Omega_h2"][1]
    target=.120
    assert 0.<ext<target<total
    q=np.asarray(s["lab_momentum_edges_GeV"])
    T=s["T_visible_GeV"]
    P=q[-1]
    theta=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    N,limits=moment_bounds(q,theta)
    assert abs(N*s["scale_factor"]**3/s["represented_species_count_a3"][1]-1.)<1e-11
    retained=total-ext
    C=coefficient(P)
    unit=N*C*P**2/(384*np.pi*m*m)
    ratios=[limits[-1]/(2*P*P),limits[-1]/(P*P)]
    stats=[]
    for gam in (1.,.98):
        H=np.sqrt(np.pi**2/45)/MPL*(T*gam)**2
        tau=unit/H
        for tag,r in zip(("equal charge populations","asymmetric charges, conservative"),ratios):
            # With arbitrary charge fractions the loss-only moment term
            # can be at most twice the charge-symmetric upper bound.
            req=mu_required((target-ext)/retained,r)
            assert abs((ext+retained*loss_floor(req,r))/target-1.)<2e-13
            t=-np.log(r)/r
            yy=solve_ivp(lambda _,x:-np.minimum(x*x,r*x),(0.,t*2+1),[1.],
                         t_eval=[t/2,t,t*2],rtol=8e-12,atol=1e-13)
            chk=max(abs(yy.y[0]/[loss_floor(x,r) for x in yy.t]-1.))
            assert chk<3e-8
            lim=[]
            for F in (0.,.001,1.,2.,40.):
                mu=tau*(1+F)**2
                p=loss_floor(mu,r)
                assert p>=max(np.exp(-r*mu),1./(1.+mu))-2e-13
                lim.append(dict(scalar_occupation_cap=F,retained_number_floor=p,
                    inclusive_fixed_count_Omega_floor=ext+retained*p,
                    same_external_support_only_floor=ext+retained/(1+mu)))
            stats.append(dict(charge_scope=tag,photon_floor=gam,support_exposure_at_zero_cap=tau,
                moment_exposure_at_zero_cap=tau*r,ratio_moment_to_support=r,transition_mu=t,
                occupation_cap_necessary_for_target=float(np.sqrt(req/tau)-1.),
                floor_ODE_independent_integration_defect=chk,cases=lim))
    return dict(input_sha256=[sha(shape),sha(audit)],temperature_GeV=T,retained_number_GeV3_at_anchor=N,
                stable_retained_p2_GeV2_lower_reconstruction_upper=limits,
                maximum_event_angle_average_pCM2_over_moment_bound=angle_identity(P),
                sampled_amp_coefficient_fraction=test(P,C),
                anchor_inclusive_omega=total,inert_exterior_omega=ext,target_omega=target,controls=stats,
                scope="Post-anchor isotropic uncorrelated forward pair annihilation only on freely redshifting retained particles. Allows selective destruction (no renormalized-moment assumption), no inverse, redistribution, gains, re-entry or initial-tail gravity correction. Fixed photon floor and a cap on future scalar occupations are assumptions. Not a bound for the general signed coupled model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shape",type=Path)
    p.add_argument("audit",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.shape,a.audit),indent=2,allow_nan=False))
