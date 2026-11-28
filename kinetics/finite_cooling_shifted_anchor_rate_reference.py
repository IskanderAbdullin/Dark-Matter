"""Local signed pair rates using the re-evaluated prompt-plus-recoil scalar."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import gross, points
from archived_collisionless_streaming_reference import Background
from finite_cooling_matched_history_los_reference import sha
from finite_cooling_signed_pair_anchor_reference import table
from finite_cooling_post_anchor_pair_diagnostic_reference import kernel
from three_species_decay_onsager_reference import m, mp


def run(eos,root):
    bg=Background(eos)
    qscale=float(bg.scale(.001))/bg.T0
    catalogs=[root/f"{part}.json" for part in ("prompt","mix")]
    records={}
    for path in catalogs:
        p=json.loads(path.read_text())
        assert p["eos_sha256"]==sha(eos)
        records[path.stem]=p["components"]
    def t(branch,case,sp,sigma,mass,g):
        v=next(x for x in records[branch] if (x["case"],x["species"],x["sigma_ln_q"])==(case,sp,sigma))
        return table(root/branch,v,qscale,mass,g)
    out=[]
    for label in ("instant","ratio5"):
        for sigma in (.1,.14):
            mxq,fs,ms=t("mix",label,"scalar",sigma,mp,1)
            q,fc,mc=t("prompt",label,"stable",sigma,m,4)
            inverse=[]
            for bins in (47,61):
                edges=np.r_[0.,np.geomspace(.000005,1.001*max(mxq,q)/qscale,bins)]
                p,_,nf=points(edges,m,fc,3)
                w=kernel(p,nf/nf.sum(),10)*(2*nf.sum())**2
                n0,u,_=gross(edges,m,mp,fc,fs,False,3,10,14,10,split_final=False)
                _,v,_=gross(edges,mp,m,fs,fc,True,3,10,14,10,split_final=False)
                assert abs(n0/(4*nf.sum())-1.)<1e-13
                assert abs(u[0]/w-1.)<2e-4
                assert u[1]/u[0]-1.<=2*ms+ms*ms+1e-12
                inverse.append(v[1]/u[1])
                out.append(dict(history=label,sigma=sigma,bins=bins,bose_dilution_ceiling=2*ms+ms*ms,
                                max_scalar=ms,max_stable=mc,quantum_vs_forward=u[1]/u[0]-1.,
                                inverse_over_forward=v[1]/u[1],unblocked_stable_rate=w))
            assert abs(inverse[0]/inverse[1]-1.)<.01
    return dict(data=out,table_catalog_sha256=[sha(p) for p in catalogs],
        scope="Frozen 1 MeV signed-rate controls at two reconstructed widths on the shifted scalar mixture and prompt stable shape. No history, continuum inverse bound, or nonlinear calculation.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("root",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.root),indent=2,allow_nan=False))
