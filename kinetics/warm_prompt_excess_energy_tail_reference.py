"""Measure the convex excess-energy effect of the resolved prompt channels."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import angle_panels, snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import Es, M, P, m, mp, width
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def tail_change(l,s,cut):
    return np.minimum(np.minimum(l,s), np.minimum(cut,np.maximum(l+s-cut,0.)))


def integrate(S,cut,n,closed,split_result=None,*,evaluators=None):
    # Optional exact phase functions are for frozen shape controls only.
    f,ep,g,el=snapshot_functions(S) if evaluators is None else evaluators
    final=ep[-1]+el[-1]
    if final<=max(cut,M):
        if split_result is not None: split_result[:]=0.
        return np.zeros(6)
    lo=max(cut,M)
    A=np.unique(np.r_[ep,cut])
    B=np.unique(np.r_[el,cut])
    def gate(wall,mass,rest):
        q=np.sqrt((wall-mass)*(wall+mass))
        hi=wall*rest+q*P
        return np.r_[M*(rest*rest+q*q)/hi,M*hi/mass**2]
    edges=np.unique(np.clip(np.r_[lo,final,gate(A,mp,Es),gate(B,m,M-Es),
                                  (A[:,None]+B).ravel()],lo,final))
    dx=np.diff(edges)/2
    z,w=roots_legendre(n)
    E=((edges[:-1,None]+edges[1:,None])/2+dx[:,None]*z).ravel()
    dE=(dx[:,None]*w).ravel()
    total=np.zeros(6)
    if split_result is not None: split_result[:]=0.
    for i in range(0,len(E),120):
        v=E[i:i+120]
        phi,W=angle_panels(v,A,B,angle=n)
        chi=v[:,None,None]-phi
        if closed:
            W=W*((phi>=ep[0])&(phi<=ep[-1])&(chi>=el[0])&(chi<=el[-1]))
        b,c=f(phi),g(chi)
        F=W*b*c
        L=W*(1+b)*(1-c)
        norm=(F+L).sum(axis=(1,2))
        h=np.divide(F.sum(axis=(1,2)),norm,out=np.zeros_like(norm),where=norm>0.)
        pos=(1.-h[:,None,None])*F
        neg=h[:,None,None]*L
        pod=pos.sum(axis=(1,2))
        np.testing.assert_allclose(pod,neg.sum(axis=(1,2)),rtol=3e-14,atol=1e-250)
        wdelta=tail_change(chi,phi,cut)
        # A binary merger raises the total tail functional; the exact
        # reverse event lowers it by the same amount.
        z=dE[i:i+120]*2/np.pi**2*np.sqrt((v-M)*(v+M))*width*M
        daughters=np.maximum(chi-cut,0.)+np.maximum(phi-cut,0.)
        if split_result is not None:
            for k,x in enumerate((phi,chi)):
                w=np.maximum(x-cut,0.)
                split_result[k]+=[z@(pos*w).sum(axis=(1,2)),z@(neg*w).sum(axis=(1,2))]
        # At a stationary parent, its excess energy cancels at each
        # fixed E_H. The independent daughter ledger avoids subtracting
        # the two large gross parent terms when both daughters are soft.
        total += [z@(v*pod),z@(pos*wdelta).sum(axis=(1,2)),z@(neg*wdelta).sum(axis=(1,2)),
                  z@((pos-neg)*wdelta).sum(axis=(1,2)),z@(pos*daughters).sum(axis=(1,2)),
                  z@(neg*daughters).sum(axis=(1,2))]
    assert 0.<=total[1]<=total[0]/2*(1+2e-14)
    tolerance=5e-14*max(*total[[1,2,4,5]],1e-240)
    assert abs((total[1]-total[2])-total[3])<=tolerance
    assert abs((total[5]-total[4])-total[3])<=tolerance
    if cut>=max(ep[-1],el[-1]):
        assert total[4]==0.
        assert not closed or total[5]==0.
    if split_result is not None:
        np.testing.assert_allclose(split_result.sum(axis=0),total[4:],rtol=5e-14,atol=1e-245)
    return total


def run(eos,press,paths):
    X=read(press)
    high=np.hypot(M,X["heavy_momentum_cut_at_late_anchor_GeV"])
    synthetic=[]
    for q in X["controls"]:
        l,s=q["light_and_scalar_energy_centers_GeV"]
        wl,ws=q["light_and_scalar_widths_GeV"]
        assert l-wl/2>high>s+ws/2
        Q=q["positive_energy_input_over_GammaM"]
        synthetic.append(dict(parent_center=q["parent_energy_center_GeV"],
                              excess_cost_over_gross_parent_energy=s/(s+l),
                              positive_excess_source_over_GammaM=Q*s/(s+l)))
    z=np.geomspace(1e-5,2.,101)
    l,s=np.meshgrid(z,z)
    for cut in (.1,.3,high):
        t=np.maximum(l+s-cut,0.)-np.maximum(l-cut,0.)-np.maximum(s-cut,0.)
        np.testing.assert_allclose(tail_change(l,s,cut),t,rtol=0.,atol=8e-16)
    out=[]
    for path in paths:
        x=read(path)
        S=x["snapshot"]
        assert x["snapshot_sha256"]==digest(S)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        T=TabulatedRadiation(eos,S["T_visible_GeV"])
        R=get_mom(p,[S["entropy_variable_"+b] for b in ("scalar","fermion")])[2].sum()
        ext=import_restart(S["positive_exterior_restart"]).energy(1.)
        H=np.sqrt((T.state(1.)[1]+R+ext)/3)/MPL
        rows=[]
        for cut in (.10,.30,.50,high):
            for closed in (True,False):
                C=np.array([integrate(S,cut,n,closed) for n in (6,10,14)])/(R*H)
                rows.append(dict(E_threshold_GeV=cut,force_all_daughters_inside=closed,
                                 orders=[6,10,14],normalized_values=C.tolist()))
        out.append(dict(snapshot_sha256=x["snapshot_sha256"],rows=rows))
    return dict(scope="Fixed energy threshold; frozen represented daughters and quasi-static prompt parent. "
                      "The positive inverse contribution is an upper on the prompt change at the same "
                      "state, not a bound over time or for other operators.",
                columns=["inverse_parent_energy_input_over_rhoH","positive_excess_energy_growth_over_rhoH",
                         "reverse_excess_energy_debit_over_rhoH","signed_gross_difference_over_rhoH",
                         "direct_daughter_excess_consumed_over_rhoH","direct_daughter_excess_emitted_over_rhoH"],
                analytic_packet_controls=synthetic,archives=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("pressure_test",type=Path)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pressure_test,a.archives),indent=2,allow_nan=False))
