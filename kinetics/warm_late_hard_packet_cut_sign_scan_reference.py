"""Scan how the initial cut current depends on the unresolved hard packet center."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import excess,split_at_cut
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,prior,archive):
    A=read(archive);S=A["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    R=TabulatedRadiation(eos,S["T_visible_GeV"])
    C=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    centers=np.unique(np.r_[np.arange(65,83)/100.,np.arange(704,723,2)/1000.])
    rows=[]
    for w in (.005,.0025):
        o=[]
        for p in centers:
            Z,*_=reactivate(S,float(p),w,cap=1.05)
            Z=split_at_cut(Z,C)
            E=get_mom(np.asarray(Z["lab_momentum_edges_GeV"]),
                       [Z["entropy_variable_"+j] for j in ("scalar","fermion")])[2].sum()
            H=np.sqrt((R.state(1.)[1]+E)/3)/MPL
            F=integrate(Z,C,10,True)/(E*H)
            W=sum(excess(Z,C))/E
            if W==0.:
                assert F[4]==0. and F[5]>=0.
            o.append(dict(hard_center_GeV=float(p),signed_initial_flux=float(F[3]),tail_excess_over_initial_energy=W,
                          consumed_and_emitted_daughter_excess=F[4:].tolist(),state_sha256=digest(Z)))
        br=[[a["hard_center_GeV"],b["hard_center_GeV"]] for a,b in zip(o[:-1],o[1:],strict=True)
            if a["signed_initial_flux"]*b["signed_initial_flux"]<0.]
        rows.append(dict(momentum_shell_halfwidth=w,sampled_sign_brackets=br,states=o))
    return dict(scope="Frozen phase choices at one late anchor; no dynamics, universal sign transition or constraint on true exterior.",
                source_sha256=A["snapshot_sha256"],energy_cut_GeV=C,cut_momenta_GeV=[float(np.sqrt(C*C-d*d)) for d in (mp,m)],
                controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive),indent=2,allow_nan=False))
