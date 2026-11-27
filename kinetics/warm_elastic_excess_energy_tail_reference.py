"""Frozen convex energy tail of the retained below-pole elastic operator."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from compact_elastic_panel_reference import compact_elastic, compile_panels
from frw_exterior_restart_reference import import_restart
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


class Panels:
    def __init__(self, state, cuts, centered, t, eta):
        self.state=state
        self.cuts=cuts
        self.values=np.zeros((len(cuts),3))
        self.n=0
        self.max_compare=0.
        self.centered=centered
        self.t,self.eta=t,eta

    def append(self,z):
        p,f,L,feature,K=z
        n=len(self.state)//2
        x=feature[:,:n]@self.state[:n]
        y=feature[:,n:]@self.state[n:]
        assert np.min(x)>0.
        block=expit(y)/(-np.expm1(-x))
        log=-x-y
        diff=log[:,None]-log[None,:]
        coeff=K*block[:,None]*block[None,:]
        signed_r=np.sign(diff)*np.exp(np.maximum(log[:,None],log[None,:]))*(-np.expm1(-abs(diff)))
        ordered=coeff*np.exp(log[:,None])
        both=coeff*signed_r
        if self.n%71==0:
            independent=feature.T@(-both.sum(axis=1))
            direct=compact_elastic(self.centered,self.t,self.eta,1.,
                                  compile_panels([z],n//2),len(self.state),True)[0]
            scale=max(np.linalg.norm(feature.T@ordered.sum(axis=1),1),1e-240)
            self.max_compare=max(self.max_compare,np.linalg.norm(independent-direct,1)/scale)
        for i,c in enumerate(self.cuts):
            # Support is closed here; a daughter tail outside the
            # retained bins is not represented by this model.
            W=np.maximum(p-c,0.)+np.maximum(f-c,0.)
            dw=W[None,:]-W[:,None]
            self.values[i]+=[np.sum(ordered*np.maximum(dw,0.)),
                             np.sum(ordered*np.maximum(-dw,0.)),
                             .5*np.sum(both*dw)]
        self.n+=1


def run(eos,pressure,archive):
    q=read(archive)
    S=q["snapshot"]
    assert q["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
    T=S["T_visible_GeV"]
    rho=get_mom(p,th)[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho+import_restart(
        S["positive_exterior_restart"]).energy(1.))/3)/MPL
    hi=np.hypot(M,read(pressure)["heavy_momentum_cut_at_late_anchor_GeV"])
    cuts=(.10,.30,.50,hi)
    cell=Collision(p)
    z=pack(th)
    model=Late(cell,H,[],True)
    qref,_,T0=model.reference(z,T)
    rows=[]
    for order in ((3,2,2,8),(4,3,3,12),(5,4,4,16)):
        tape=Panels(z,cuts,z-qref,T0,model.eta)
        setup(max(T0,T),cell,T0,model.eta,*order,32.,.82,cache=tape,cache_only=True,
              energy_observable_cuts=cuts)
        vals=tape.values/(rho*H)
        scale=np.maximum(1e-245,np.max(abs(vals[:,:2]),axis=1))
        assert np.max(abs(vals[:,0]-vals[:,1]-vals[:,2])/scale)<2e-12
        assert np.all(vals[np.asarray(cuts)>=max(np.hypot(k.mass,p[-1]) for k in cell.c)]==0.)
        assert tape.max_compare<3e-8
        rows.append(dict(order=list(order),lab_panels=tape.n,rate_over_retained_rhoH=vals.tolist(),
                         max_projected_panel_error_over_gross_against_compact=tape.max_compare))
    return dict(scope="Frozen retained-only elastic channel with strict subpole and compact model lab support. "
                      "No exterior elastic output or higher pole channels. Not an integral or upper bound.",
                snapshot_sha256=q["snapshot_sha256"],cuts_GeV=cuts,
                columns=["positive_gross_excess","negative_gross_excess","signed_excess"],controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","pressure","archive"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pressure,a.archive),indent=2,allow_nan=False))
