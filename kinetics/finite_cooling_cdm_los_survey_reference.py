"""CDM-only metric driver setup and low-memory fixed-metric tracer summary."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from finite_cooling_pair_split_tracer_reference import fields

K=np.array([1.,2.,4.,8.,10.,16.,22.])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(seed,out,modes=None):
    out.mkdir(parents=True,exist_ok=True)
    modes=K if modes is None else np.asarray(modes,dtype=float)
    assert len(modes)>0 and np.all(np.diff(modes)>0) and .5<=min(modes)<=max(modes)<=60
    L=fields(seed)
    assert L["gauge"]=="synchronous" and abs(float(L["h"])-.6736)<1e-12
    d={k:v for k,v in L.items() if "ncdm" not in k and k!="evolver"}
    d.update({"root":str(out/"cdm_los_"),"z_pk":"0,5,20","P_k_max_h/Mpc":f"{max(modes):g}",
        "k_min_tau0":"2500","k_per_decade_for_pk":"6","write background":"no","perturbations_verbose":"3",
        "k_output_values":", ".join(str(.6736*k) for k in modes)})
    # Keep general formatting stable and ensure no ncdm fields at all:
    assert all("ncdm" not in k for k in d)
    cfg=out/"cdm_los.ini"
    cfg.write_text("".join(f"{k} = {v}\n" for k,v in d.items()))
    return dict(seed_ini_sha256=sha(seed),k_h_per_Mpc=modes.tolist(),native_parameters={k:v for k,v in d.items() if k!="root"},
        scope="CDM+baryon+radiation metric with no scalar feedback. The synthetic tracer enters only in the separately evaluated LOS, not in this native CLASS background.")


def summarize(records,reference,fluid,split):
    base=json.loads(reference.read_text())
    cmp={round(x["k_h_per_Mpc"],10):x["cases"][-1]["los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"][0][0] for x in base["results"]}
    F=json.loads(fluid.read_text())
    S=json.loads(split.read_text())
    rkos=next(x["summed_tracer_delta_over_cdm"] for x in S["cases"] if x["case"]=="split_off_l36_rk")
    lastK=[]
    rows=[]
    for path in records:
        A=json.loads(path.read_text())
        z=float(A["z"])
        lastK=[float(x["k_h_per_Mpc"]) for x in A["results"]]
        assert np.max(abs(np.array(lastK)-K))<1e-12
        key="los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"
        vals=np.array([x["cases"][-1][key] for x in A["results"]])
        C=vals[:,:,0]
        rho=np.array([r["cases"][-1]["summed_prompt_recoil_with_full_IC"] for r in A["results"]])
        assert np.max(abs(rho-C[:,0]))<3e-7
        w=next(x for x in F["data"] if x["z"]==z)
        kk=w["probes"].index(22.)
        row=dict(z=z,k_h_per_Mpc=K.tolist(),joined_prompt_recoil=C.T.tolist(),
                largest_refinement_difference=float(max(x["max_time_and_momentum_refinement_difference"] for x in A["results"])),
                tested_joined_samples_below_half=int(np.count_nonzero(C[:,0]<.5)),
                max_joined_vs_energy_weighted_split=float(np.max(abs(C[:,0]-rho))),
                default_fluid_comparator_k22=w["wide_higher_delta_over_cdm"][kk])
        if z==20.:
            ix=np.searchsorted(K,[4.,10.])
            d=[float(y-cmp[float(k)]) for y,k in zip(C[ix,0],K[ix])]
            s=(C[ix,0]-rkos).tolist()
            assert max(abs(np.array(d)))<1e-6
            row.update(anchor_k_h_per_Mpc=[4.,10.],absolute_cdm_metric_minus_native_recorded_metric=d,
                       absolute_cdm_metric_los_minus_native_split_rk=s)
        rows.append(row)
    assert sorted(r["z"] for r in rows)==[0.,5.,20.]
    return dict(data=rows,cdm_metric_los_records_sha256=[sha(p) for p in records],reference_sha256=[sha(p) for p in (reference,fluid,split)],
        scope="Discrete LOS tracer samples over a CDM-only native fixed metric, not a self-consistent mixture transfer, observed small-scale phenomenon or continuously certified threshold.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--seed",type=Path)
    p.add_argument("--modes",type=float,nargs="+")
    p.add_argument("--out",type=Path)
    p.add_argument("--reference",type=Path)
    p.add_argument("--fluid",type=Path)
    p.add_argument("--split",type=Path)
    p.add_argument("records",type=Path,nargs="*")
    a=p.parse_args()
    result=prepare(a.seed,a.out,a.modes) if a.seed is not None else summarize(a.records,a.reference,a.fluid,a.split)
    print(json.dumps(result,indent=2,allow_nan=False))
