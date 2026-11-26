"""Finite-volume Pauli-bounded χℓ 2->2 NONLINEAR time pilot (synthetic).

Recompute exact same 5D on-shell finite-vector Pauli event source from
current nodal occupations, reconstructed via consistent Fermi logits.
Node particle capacity = exact radial density-of-states cell integral:
 w_i = g/(6*pi²) [(E_hi²-m²)^(3/2)-(E_lo²-m²)^(3/2)].
Then df_i/dt = dN_i/w_i.

Adaptive Euler step-doubling: accept TWO half steps only when EVERY
intermediate f lies strictly (0,1), entropy is nondecreasing, and full
three-quantity (all chi/lepton occupations) local difference <= atol.
Rejected proposals halve dt; NO Fermi clipping or thermal reset.
Preserves each χ/ℓ number and total discrete energy to floating error.
Step-doubling estimator is NOT a rigorous GLOBAL time-error bound;
accepted method remains first order, not a physical χ gamma_el.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit, xlogy
from kinetics.paired_chi_lepton_barycentric_5d import projected_pair
from kinetics.paired_chi_lepton_logit_galerkin_5d import logit_interpolant


def phase_space_cell_weights(nodes, mass, degeneracy=4):
    """Physical radial spin-weighted d³p cells, clipped at Emax and mass."""
    nodes=np.asarray(nodes,float)
    if (mass<0 or degeneracy<=0 or len(nodes)<3 or
        not np.all(np.isfinite(nodes)) or
        not np.all(np.diff(nodes)>0) or
        abs(nodes[0]-mass)>1e-12 or
        not np.allclose(np.diff(nodes),nodes[1]-nodes[0],rtol=1e-12)):
        raise ValueError("finite uniform energy grid starting at species rest mass")
    h=nodes[1]-nodes[0]
    lo=np.maximum(mass,nodes-h/2)
    hi=np.minimum(nodes[-1],nodes+h/2)
    w=degeneracy/(6*math.pi**2)*(
        np.maximum(hi*hi-mass*mass,0)**1.5 -
        np.maximum(lo*lo-mass*mass,0)**1.5)
    if np.any(w<=0):
        raise ArithmeticError("zero phase-space capacity")
    return w


def quantum_entropy(f,w):
    """Fermi entropy including well-defined exact 0/1 limits."""
    if (np.any(f<0) or np.any(f>1) or
        not np.all(np.isfinite(f))):
        raise ValueError("Fermi-bound violation")
    return -float(np.dot(w,xlogy(f,f)+xlogy(1-f,1-f)))


def system(M=.1,m=0.,Te=1.,cut=12.,bins=10,mode="oscillatory"):
    if not (M>0 and m>=0 and Te>0 and cut>0
            and isinstance(bins,int) and bins>=4):
        raise ValueError
    hi=M+m+2*cut*Te
    cg=np.linspace(M,hi,bins+1)
    lg=np.linspace(m,hi,bins+1)
    if mode=="oscillatory":
        fc=expit(-cg+.3*np.sin(4*cg))
        fl=expit(-lg+.3*np.cos(4*lg))
    elif mode=="smooth":
        fc=expit(-cg/.9+.15*np.exp(-cg/.9))
        fl=expit(-lg+.05*np.exp(-lg))
    elif mode=="equal_temperature":
        fc=expit(.12-cg)
        fl=expit(-.08-lg)
    else:
        raise ValueError("unknown initial mode")
    return {"M":M,"m":m,"Te":Te,"cut":cut,"cg":cg,"lg":lg,
            "wc":phase_space_cell_weights(cg,M,4),
            "wl":phase_space_cell_weights(lg,m,4),
            "initial_chi":fc,"initial_lepton":fl,
            "historical_gchi_glepton_4each_not_physical_nu_flavors":True}


def ledger(cfg,fc,fl):
    return {"Nchi":float(np.dot(cfg["wc"],fc)),
            "Nlepton":float(np.dot(cfg["wl"],fl)),
            "total_energy":float(np.dot(cfg["wc"],cfg["cg"]*fc)+
                                 np.dot(cfg["wl"],cfg["lg"]*fl)),
            "entropy":quantum_entropy(fc,cfg["wc"])+
                      quantum_entropy(fl,cfg["wl"])}


def rhs(cfg,fc,fl):
    if (np.any(fc<=0) or np.any(fc>=1) or np.any(fl<=0) or
        np.any(fl>=1)):
        raise ValueError("strict Fermi interior for logit reconstruction")
    f=logit_interpolant(cfg["cg"],fc)
    g=logit_interpolant(cfg["lg"],fl)
    r=projected_pair(f,g,M=cfg["M"],m=cfg["m"],Te=cfg["Te"],
                     cut_kinetic_T=cfg["cut"],bins=len(fc)-1)
    return (r["dNchi_bins"]/cfg["wc"],
            r["dNlepton_bins"]/cfg["wl"],
            r["direct_event_entropy"],
            r["entropy_identity_abs_residual"]
                if "entropy_identity_abs_residual" in r else
                abs(r["direct_event_entropy"]-r["interpolated_bin_entropy"]))


def euler_if_physical(cfg,fc,fl,dt):
    dc,dl,sdot,identity_error=rhs(cfg,fc,fl)
    xc=fc+dt*dc
    xl=fl+dt*dl
    if (np.any(xc<=0) or np.any(xc>=1) or np.any(xl<=0) or
        np.any(xl>=1)):
        return None
    before=ledger(cfg,fc,fl)
    after=ledger(cfg,xc,xl)
    # No entropy post-processing or resetting of the occupations.
    if after["entropy"]<before["entropy"]-1e-12:
        return None
    return xc,xl,sdot,identity_error


def evolve(*,bins=10,mode="oscillatory",total_time=1.e10,
           initial_dt=1.e10,atol=1.e-7,max_attempts=300):
    if not (total_time>0 and initial_dt>0 and atol>0 and
            all(math.isfinite(x) for x in (total_time,initial_dt,atol))):
        raise ValueError
    cfg=system(bins=bins,mode=mode)
    fc=cfg["initial_chi"].copy()
    fl=cfg["initial_lepton"].copy()
    first=ledger(cfg,fc,fl)
    t=0.;dt=initial_dt
    accepted=rejected=0
    errors=[]
    entropy_increments=[]
    accepted_dt=[]
    for attempt in range(max_attempts):
        if t>=total_time*(1-1e-14):
            break
        dt=min(dt,total_time-t)
        full=euler_if_physical(cfg,fc,fl,dt)
        half=euler_if_physical(cfg,fc,fl,dt/2)
        two=None if half is None else euler_if_physical(
            cfg,half[0],half[1],dt/2)
        if full is not None and two is not None:
            err=max(float(np.max(abs(full[0]-two[0]))),
                    float(np.max(abs(full[1]-two[1]))))
            if err<=atol:
                before=ledger(cfg,fc,fl)
                fc,fl=two[0].copy(),two[1].copy()
                after=ledger(cfg,fc,fl)
                entropy_increments.append(after["entropy"]-before["entropy"])
                accepted_dt.append(dt)
                errors.append(err)
                accepted+=1;t+=dt
                dt*=min(2.,max(1.1,.9*math.sqrt(atol/max(err,1e-30))))
                continue
        rejected+=1
        dt*=.5
        if dt<1e-18:
            raise ArithmeticError("positivity/entropy error test cannot find physical step")
    else:
        raise RuntimeError("maximum adaptive attempts exhausted")
    final=ledger(cfg,fc,fl)
    return {"mode":mode,"bins":bins,
        "numerical_time_not_cosmological":total_time,
        "accepted_steps":accepted,"rejected_proposals":rejected,
        "min_accepted_dt":min(accepted_dt),"max_accepted_dt":max(accepted_dt),
        "max_step_doubling_spectral_abs_error":max(errors),
        "number_chi_error":final["Nchi"]-first["Nchi"],
        "number_lepton_error":final["Nlepton"]-first["Nlepton"],
        "total_energy_error":final["total_energy"]-first["total_energy"],
        "entropy_gain":final["entropy"]-first["entropy"],
        "min_accepted_step_entropy_gain":min(entropy_increments),
        "min_fchi":float(min(fc)),"max_fchi":float(max(fc)),
        "min_flepton":float(min(fl)),"max_flepton":float(max(fl)),
        "fchi_first_bin":float(fc[0]),"flepton_first_bin":float(fl[0]),
        "adaptive_step_doubling_not_rigorous_global_error_bound":True,
        "first_order_Euler_two_half_steps":True,
        "class_ready":False}


if __name__=="__main__":
    import json
    results={}
    for mode in ("oscillatory","smooth"):
        results[mode]={str(a):evolve(mode=mode,atol=a)
                        for a in (1e-6,1e-7,1e-8)}
    print(json.dumps(results,indent=2))
