"""Matched difference of relative time shifts between q16 and q32."""
import argparse
import json
import gzip
from pathlib import Path

import numpy as np

from warm_early_serial_support_audit_reference import rows, sha


def read(p):
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix==".gz" else p.read_bytes())


def run(archive, seed16, seed32, a16, a32, b16, b32, prefix=False,base_blocks=4):
    assert base_blocks in (4,8)
    files=(archive,seed16,seed32,a16,a32,b16,b32)
    A=read(archive)["cases"]
    S=[read(p)["snapshots"][-1] for p in (seed16,seed32)]
    s=[sha(p) for p in (seed16,seed32)]
    P=[read(p) for p in (a16,a32,b16,b32)]
    assert abs(S[0]["T_visible_GeV"]/S[1]["T_visible_GeV"]-1.)<1e-10
    for i,B in enumerate(P):
        grid=(base_blocks,2*base_blocks)[i//2]
        assert B["blocks_total"]==grid
        if i!=3 or not prefix:
            assert B["blocks_finished"]==grid
        for key in ("schema","eos_sha256","local_prompt_jacobian","reference_sha256"):
            assert B["plan"][key]==P[0]["plan"][key]
    for j in range(2):
        assert np.max(abs(np.asarray(P[j+2]["plan"]["outer_log_scale_grid"])[::2]-
                          P[j]["plan"]["outer_log_scale_grid"]))<1e-15
    R=[rows(S[i%2],p,A[i%2],s[i%2]) for i,p in enumerate(P)]
    out=[]
    for h in zip(R[0],R[1],R[2][1::2],R[3][1::2]):
        T=h[0]["T_GeV"]
        b=h[0]["b"]
        assert max(abs(x["b"]/b-1.) for x in h)<1e-13
        assert max(abs(x["T_GeV"]/T-1.) for x in h)<1e-12
        V=np.asarray([x["inclusive_fixed_count_omega"] for x in h])
        L=np.log(V).reshape(2,2,2)
        support=L[:,1]-L[:,0]
        time=L[1]-L[0]
        D=np.diff(support,axis=0)[0]
        check=np.max(abs(D-(time[1]-time[0])))
        assert check<2e-14
        out.append(dict(b=b,T_GeV=T,omegas_4q16_4q32_8q16_8q32=V.tolist(),
                        log_support_shift_4_and_8=support.tolist(),
                        log_time_shift_q16_and_q32=time.tolist(),
                        time_support_interaction=D.tolist(),interaction_to_four_block_support=abs(D/support[0]).tolist(),
                        support_sign_agrees=(support[0]*support[1]>0).tolist()))
    if base_blocks!=4:
        for x in out:
            for source,dest in (("omegas_4q16_4q32_8q16_8q32","omegas_8q16_8q32_16q16_16q32"),
                                ("log_support_shift_4_and_8","log_support_shift_8_and_16"),
                                ("interaction_to_four_block_support","interaction_to_eight_block_support")):
                x[dest]=x.pop(source)
    result=dict(complete=len(out)==base_blocks,matched=out,input_sha256=[sha(p) for p in files],
         scope="Difference-in-differences of terminal inclusive rest-count proxies on exactly matched q16 and q32 time grids. Same sign across tested resolutions is not continuum convergence or a general physical error bar.")
    if base_blocks!=4:
        result["coarse_and_refined_steps"]=[base_blocks,2*base_blocks]
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("archive","seed16","seed32","a16","a32","b16","b32"):
        p.add_argument(key,type=Path)
    p.add_argument("--allow-prefix",action="store_true")
    p.add_argument("--base-blocks",type=int,default=4)
    a=p.parse_args()
    print(json.dumps(run(a.archive,a.seed16,a.seed32,a.a16,a.a32,a.b16,a.b32,a.allow_prefix,a.base_blocks),indent=2,allow_nan=False))
