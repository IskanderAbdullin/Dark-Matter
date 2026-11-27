"""Check the local sign boundary against hard packet shell width."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import excess,split_at_cut
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,prior,archive):
    A=read(archive)
    S=A["snapshot"]; assert digest(S)==A["snapshot_sha256"]
    T=TabulatedRadiation(eos,S["T_visible_GeV"])
    C=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    momentum_cut=np.sqrt((C-m)*(C+m))

    def kernel(p,w,n):
        D,*_=reactivate(S,p,w,cap=1.05)
        D=split_at_cut(D,C)
        N,K,E=get_mom(np.asarray(D["lab_momentum_edges_GeV"]),
                       [D["entropy_variable_"+q] for q in ("scalar","fermion")])
        rho=E.sum()
        H=np.sqrt((T.state(1.)[1]+rho)/3)/MPL
        z=integrate(D,C,n,True)/(rho*H)
        return float(z[5]-z[4]),float(sum(excess(D,C))/rho)

    rows=[]
    for w in (.005,.0025,.001,.00025):
        first,last=.708,.718
        f0=kernel(first,w,8)[0]; f1=kernel(last,w,8)[0]
        assert f0>0.>f1
        root=brentq(lambda p:kernel(p,w,8)[0],first,last,xtol=3e-10)
        at=[kernel(root,w,n) for n in (8,12,16)]
        dp=2e-6
        J=[kernel(root+q*dp,w,12)[0] for q in (-1,1)]
        assert J[0]>0.>J[1]
        rows.append(dict(physical_hard_shell_halfwidth_GeV=w,local_crossing_candidate_GeV=float(root),
                         sampled_sign_bracket_GeV=[float(root-dp),float(root+dp)],
                         initial_signed_flux_at_bracket=J,
                         signed_flux_and_excess_at_root_by_order=[list(v) for v in at],
                         shift_below_light_cut_in_halfwidths=float((momentum_cut-root)/w)))
    return dict(scope="Local sign changes in one family of projected frozen packet shapes only. "
                      "Not a boundary for arbitrary exterior distributions or a duration estimate.",
                source_sha256=A["snapshot_sha256"],cut_energy_GeV=C,
                light_cut_momentum_GeV=momentum_cut,controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive),indent=2,allow_nan=False))
