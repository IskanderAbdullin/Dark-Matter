"""Audit the sampled near-cut histories without interpolating sign changes."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(*folders):
    outputs=[]
    anchor=None
    Ns=[]
    hs=[]
    eroot=[]
    paired_moments=[]
    for folder in folders:
        z=read(folder/"summary.json")
        p=z["plan"]
        attrs=(p["archive_sha256"],p["eos_sha256"],p["normalized_cut_GeV"],
               tuple(x["time"] for x in z["rows"]))
        if anchor is None:
            anchor=attrs
        assert anchor==attrs
        first=read(folder/"step-0.json")["snapshot"]
        q=np.asarray(first["lab_momentum_edges_GeV"])
        theta=[first["entropy_variable_"+n] for n in ("scalar","fermion")]
        N,_,E=get_mom(q,theta)
        Ns.append(N.sum(axis=1))
        hs.append(z["H_GeV"])
        rho=E.sum()
        eroot.append(float(rho))
        paired_moments.append((p["hard_GeV"],p["shell_halfwidth_GeV"],N.sum(axis=1),E.sum(axis=1)))
        rows=z["rows"]
        for j,x in enumerate(rows):
            r=read(folder/f"step-{j}.json")
            y=r.pop("snapshot")
            assert r==x and r["hash"]==digest(y) and r["plan"]==p
            assert abs(sum(excess(y,p["normalized_cut_GeV"]))/rho/r["tail_energy"]-1.)<5e-11
            assert max(r["discrete_N_and_E_errors"])<8e-12
            assert len(r["substep_solver"])==(0 if j==0 else p["substeps"])
        brackets=[]
        for a,b in zip(rows,rows[1:]):
            if a["signed_native_flux"]*b["signed_native_flux"]<0.:
                brackets.append([a["time"],b["time"]])
        output=dict(hard_GeV=p["hard_GeV"],halfwidth_GeV=p["shell_halfwidth_GeV"],n=p["substeps"],seed=p["seed_sha256"],
                    W=[q["tail_energy"] for q in rows],rates=[q["signed_native_flux"] for q in rows],
                    times=[q["time"] for q in rows],sampled_sign_brackets=brackets,
                    Wmax=float(max(q["tail_energy"] for q in rows)),
                    W_final_over_initial=float(rows[-1]["tail_energy"]/rows[0]["tail_energy"]))
        outputs.append(output)
    assert max(np.max(abs(N/Ns[0]-1.)) for N in Ns)<8e-13
    matched_moments=[]
    for i,(q,w,n,e) in enumerate(paired_moments):
        for Q,W,N,E in paired_moments[i+1:]:
            if Q==q and W!=w:
                defect=max(float(max(abs(n/N-1.))),float(max(abs(e/E-1.))))
                assert defect<4e-13
                matched_moments.append(dict(hard_GeV=q,halfwidths=[w,W],
                            per_species_N_and_E_bitwise_equal=bool(np.array_equal(n,N) and np.array_equal(e,E)),
                            max_relative_moment_difference=defect))
    return dict(same_center_width_moment_checks=matched_moments,scope="Exact archive and cut alignment; recorded adjacent sign changes only, no fitted dynamics.",
                rows=outputs,normalization_H_max_relative_shift=max(abs(h/hs[0]-1.) for h in hs),
                max_initial_total_daughter_energy_relative_shift=max(abs(e/eroot[0]-1.) for e in eroot))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folders",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(*a.folders),indent=2,allow_nan=False))
