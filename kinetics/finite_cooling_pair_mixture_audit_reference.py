"""CDF and geodesic bracket audit for the prompt + independent pair mixture."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from archived_collisionless_streaming_reference import Background
from finite_cooling_reheat_distribution_reference import stream_bounds
from three_species_decay_onsager_reference import mp


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(eos,prompt,pair):
    a,b=[json.loads(p.read_text()) for p in (prompt,pair)]
    esh=sha(eos)
    assert a["eos_sha256"]==b["eos_sha256"]==esh
    B=Background(eos)
    cut=np.array([.05,.1,.2,.5,1.,2.,4.,8.,16.])
    edges=np.r_[b["recoil_p_over_TF_lower_edges"],np.inf]
    ia=np.searchsorted(a["cdf_thresholds_p_over_T"],cut)
    ib=np.searchsorted(b["recoil_shape_probe_thresholds"],cut)
    assert np.array_equal(np.asarray(a["cdf_thresholds_p_over_T"])[ia],cut)
    assert np.array_equal(np.asarray(b["recoil_shape_probe_thresholds"])[ib],cut)
    ans=[]
    for p,q in zip(a["cases"],b["episodes"]):
        assert p["T_R_GeV"]==q["T_R_GeV"]
        x,y=q["modes"]
        sx,sy=[z["new_scalar_recoil_shape"] for z in (x,y)]
        S=sy
        n=np.sum(S["count_fraction_by_epoch"],axis=0)
        p2=np.sum(S["q2_weighted_count_fraction_by_epoch"],axis=0)
        assert abs(n.sum()-1.)<3e-12
        streams=stream_bounds(B,edges,n,p2,mp,8)
        old=p["species"][0]
        new=float(np.sum(np.array(y["new_scalar_omega"])[:,0]))
        N=old["Omega_h2"]+new
        first=new/N
        c0=np.asarray(old["CDF_at_probes"])[ia]
        c1=np.asarray(S["fraction_below_thresholds"])[ib]
        maxdiff=np.max(np.abs(np.array(sx["fraction_below_thresholds"])-np.array(sy["fraction_below_thresholds"])))
        assert maxdiff<.012
        mix=c0*(1-first)+c1*first
        assert np.all(np.diff(mix)>=0.) and 0<=mix.min()<=mix.max()<=1.
        bb=np.asarray(old["streaming"]["mean_Mpc"])*(1-first)+np.asarray(streams["mean_Mpc"])*first
        ans.append(dict(T_R_GeV=p["T_R_GeV"],Tmax_GeV=p["Tmax_GeV"],
           conditional_mixture_scalar_Omega_h2=N,pair_fraction_of_two_component_count=first,
           cdf_thresholds_p_over_T=cut.tolist(),prompt_CDF=c0.tolist(),pair_CDF=c1.tolist(),
           first_order_mixture_CDF=mix.tolist(),pair_shape_absolute_max_CDF_probe_difference=float(maxdiff),
           standalone_pair_recoil_streaming=streams,mixture_mean_streaming_Mpc=bb.tolist()))
    return dict(eos_sha256=esh,prompt_sha256=sha(prompt),pair_sha256=sha(pair),cases=ans,
        scope="Convex combination of represented prompt and passive first-order pre-anchor pair sources, no further reactions. These are conditional geodesic statistics, not perturbative power or coupled physical yields.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","prompt","pair"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prompt,a.pair),indent=2,allow_nan=False))
