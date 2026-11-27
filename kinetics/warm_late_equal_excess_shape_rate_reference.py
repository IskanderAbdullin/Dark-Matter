"""Test whether archive moments plus the initial convex tail determine its sign."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import excess, split_at_cut
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos, prior, source):
    data=read(source)
    S=data["snapshot"]
    assert digest(S)==data["snapshot_sha256"]
    cut=float(np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"]))
    T=TabulatedRadiation(eos,S["T_visible_GeV"])
    target=2e-16

    def selected(center, width):
        A,_,_,_,_=reactivate(S,center,width,cap=1.05)
        A=split_at_cut(A,cut)
        p=np.asarray(A["lab_momentum_edges_GeV"])
        N,K,E=get_mom(p,[A["entropy_variable_"+n] for n in ("scalar","fermion")])
        rho=E.sum()
        H=np.sqrt((T.state(1.)[1]+rho)/3.)/MPL
        W=np.asarray(excess(A,cut))
        return W.sum()/rho,A,W/rho,N.sum(axis=1),E.sum(axis=1),rho*H

    result=[]
    for width,bracket in ((.0025,(.711,.716)),(.00025,(.7145,.716))):
        a,b=bracket
        assert selected(a,width)[0]<target<selected(b,width)[0]
        center=brentq(lambda center: selected(center,width)[0]/target-1.,a,b,xtol=1e-15)
        total,A,parts,N,E,norm=selected(center,width)
        assert abs(total/target-1)<2e-11
        rates=[]
        for order in (8,12,16):
            split=np.empty((2,2))
            v=integrate(A,cut,order,True,split_result=split)/norm
            split/=norm
            assert abs(v[3]-(split[:,1]-split[:,0]).sum()) < 6e-14*split.sum()
            rates.append(dict(order=order,signed=float(v[3]),gross=float(v[1]+v[2]),
                              species_consumption_and_emission=split.tolist()))
        assert max(abs(r["signed"]-rates[-1]["signed"]) for r in rates)<2e-5*rates[-1]["gross"]
        result.append(dict(center_GeV=center,halfwidth_GeV=width,root_bracket_GeV=bracket,
                       W_over_rho=float(total),per_species_W_over_rho=parts.tolist(),
                       species_N=N.tolist(),species_E=E.tolist(), angular=rates,
                       snapshot_sha256=digest(A),snapshot=A))
    diff=max(abs(np.r_[result[0]["species_N"],result[0]["species_E"]]/
                 np.r_[result[1]["species_N"],result[1]["species_E"]]-1.))
    assert diff<1e-13
    assert all(a["signed"]>0>b["signed"] for a,b in zip(result[0]["angular"],result[1]["angular"],strict=True))
    return dict(scope="Prescribed reconstruction, matched scalar-plus-light W and bulk moments. "
                      "Frozen closed prompt current only; per-species W need not agree.",
                source_sha256=data["snapshot_sha256"],cut_GeV=cut, W_over_rho_target=target,
                max_relative_species_moment_difference=float(diff),cases=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","source"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.source),indent=2,allow_nan=False))
