"""Audit decay-mode growth signs and the minimum daughter counts needed to reverse them."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from early_pole_decay_mode_clock_reference import intervals,parent_integral
from early_pole_shell_hartree_feedback_reference import S,Q


def stats(M,m,b,p,f,edges,g=4.):
    E,lo,hi,span,rest=intervals(M,m,b,p)
    K=parent_integral(lo,hi)
    if K<1e-12:
        return None
    l=np.maximum(lo,E+edges[:-1])
    h=np.minimum(hi,E+edges[1:])
    D=np.maximum(h-l,0.)
    G=np.log1p(np.exp(-l)*(-np.expm1(-D))/(1+np.exp(-l-D)))
    v=G.sum()
    assert abs(v/K-1)<.002, (M,m,b,p,v,K)
    minE=max(lo-E,b)
    lo_sq=(minE-b)*(minE+b)
    d2=K*(2*minE+K)
    a=np.sqrt(lo_sq)
    z=np.sqrt(lo_sq+d2)
    # Counts per energy increase monotonically; fill the cheapest
    # interval to f=1 until its unweighted length reaches K.
    bound=g*(d2/(z+a))*(z*z+z*a+a*a)/(6*np.pi**2*S)
    body=f@D
    seed=(G@(1-f))/K
    return dict(scalar_q_over_T=p,fermion_min_energy_over_T=minE,parent_E_lower=lo,
                parent_E_upper=hi,mean_parent_for_sign_flip=K/span,
                measured_mean_daughter_over_mode_window=body/span,
                B_over_I=body/K,remaining_seed_fraction_of_empty=seed,
                minimum_global_daughter_Y_for_sign_flip=bound)


def control():
    # Independent direct number integral from the extremal filled band.
    M=.05;m=.0011285;p=Q
    E,lo,hi,span,rest=intervals(M,m,0.,p)
    K=quad(lambda x:expit(-x),lo,hi,epsabs=1e-13)[0]
    v=quad(lambda x:x*x,lo-E,lo-E+K,epsabs=1e-13)[0]*4/(2*np.pi**2*S)
    q=stats(M,m,0.,p,np.array([0.]),np.array([0.,30.]))
    assert abs(v/q['minimum_global_daughter_Y_for_sign_flip']-1)<3e-12
    return q


def run(root):
    out=[]
    for label in ('small','fine','small-direct','fine-direct','vacuum-direct','only-pauli','finite-daughter'):
        path=root/(label+'.json')
        if not path.exists() or not path.stat().st_size:
            continue
        u=json.loads(path.read_text())
        M=u['checkpoints'][-1]['x'];m=u['checkpoints'][-1]['scalar_pole_over_T']; b=u['fermion_pole_over_T']
        ed=np.hypot(np.array(u['fermion_p_edges']),b)
        f=np.array(u['final_daughter_per_spin_and_charge_occupations'])
        se=np.array(u['scalar_p_edges'])
        counts=np.array(u['final_scalar_n_per_visible_entropy'])
        p=np.sqrt(se[:-1]*se[1:])
        peak=float(p[np.argmax(counts/np.diff(np.log(se)))])
        y=u['checkpoints'][-1]['Y']
        data=[stats(M,m,b,p,f,ed) for p in (peak,Q,1e-5)]
        for row in [d for d in data if d and d['B_over_I']>=1.]:
            assert y>=row['minimum_global_daughter_Y_for_sign_flip']*(1-1e-5)
        out.append(dict(name=label,x=M,m=m,b=b,p_over_T_with_most_scalar_count_per_log_bin=peak,
                        final_total_daughter_Y_by_pair_conservation=y,mode_tests=data))
    return dict(scope="Instantaneous finite energy windows in the retained held-parent model. "
                      "The Pauli minimum is necessary, not sufficient; arbitrary daughter shape allowed "
                      "by covering the lowest energetic states first. Final binned occupations used "
                      "only as a piecewise constant diagnostic.",normalization_support_control=control(),tests=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path)
    print(json.dumps(run(p.parse_args().runs),indent=2,allow_nan=False))
