"""Frozen-loss count tangent for the independent low-reheat examples."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_portal_reference import SourceFold
from three_species_decay_onsager_reference import m, mp


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(eos, prompt, shape, post, feedback):
    files=[prompt,shape,post,feedback]
    P,S,A,D=[json.loads(p.read_text()) for p in files]
    assert all(d["eos_sha256"]==sha(eos) for d in (P,S,A,D))
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    target=P["target_stable_Omega_h2"]
    assert target==.120
    controls=[]
    for order,n,q in (((7,18,14),130,4),((10,26,20),230,6)):
        K=SourceFold(eos,order,n,q)
        rows=[]
        for i in range(2):
            a,s,d=A["cases"][i],S["episodes"][i],D["cases"][i]
            assert a["T_R_GeV"]==s["T_R_GeV"]==d["T_R"] and a["Tmax_GeV"]==s["Tmax_GeV"]==d["T_max"]
            ref=s["T_R_GeV"]
            ratio=s["Tmax_GeV"]/ref
            assert ratio==(1.,5.)[i]
            def stable(T):
                return f*m*(K.post(T)+K.before(T,ratio))[1]
            root=brentq(lambda T:stable(T)-target,.005,.010,xtol=2e-14)
            slopes=[]
            for epsilon in (1e-4,5e-5):
                h=epsilon*root
                slopes.append((stable(root+h)-stable(root-h))/(2*h))
            assert abs(slopes[0]/slopes[1]-1.)<3e-6
            gross=np.sum(np.asarray(s["modes"][-1]["new_scalar_omega"])[:,0])
            fb=sum(d["feedback_midpoint_converted_omega_scalar"])
            a1=a["controls"][-1]
            late=a1["perturbative_unblocked_scalar_Omega_h2_equivalent"]
            tau=a1["integrated_vacuum_forward_log_expansion_optical_depth"]
            assert abs(late/(tau*target*mp/m)-1.)<2e-7
            # The alternate early fold uses its own matched fixed nodes.
            # It is a graph-resolution/stable-loss check, not combined
            # simultaneously with the separate passive inverse estimate.
            ans=[]
            for label,pre in (("independent_forward",gross),("fixed_node_depletion_check",fb)):
                loss=(pre+late)*m/mp
                assert 0.<loss/target<.002
                corr=brentq(lambda T:stable(T)-target-loss,.005,.010,xtol=2e-14)
                actual=corr-root
                tangent=loss/slopes[-1]
                assert abs(actual/tangent-1.)<5e-4
                ans.append(dict(case=label,stable_omega_depletion=loss,early_new_scalar_rest_omega=pre,
                                late_new_scalar_rest_omega=late,shifted_root_GeV=corr,delta_T_keV=actual*1e6,
                                tangent_delta_T_keV=tangent*1e6))
            rows.append(dict(T_max_over_T_R=ratio,archived_unconverted_root_GeV=ref,recomputed_unconverted_root_GeV=root,
                             stable_count_temperature_slope_per_GeV=slopes[-1],models=ans))
        controls.append(dict(source_order=order,log_grid=n,panel_order=q,rows=rows))
    for a,b in zip(*(z["rows"] for z in controls)):
        for x,y in zip(a["models"],b["models"]):
            assert abs(x["delta_T_keV"]-y["delta_T_keV"])<.001
    return dict(controls=controls,input_sha256={p.name:sha(p) for p in [eos,*files]},
                scope="A local scalar/stable count conversion and fixed-loss finite difference on the prescribed portal clock. The early and late forward yields do not overlap at 1 MeV. This keeps both losses fixed when shifting T_R, omits a jointly iterated inverse source and feedback, and does not generate a new spectrum or infer a physical reheating temperature.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("eos","prompt","shape","post","feedback"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prompt,a.shape,a.post,a.feedback),indent=2,allow_nan=False))
