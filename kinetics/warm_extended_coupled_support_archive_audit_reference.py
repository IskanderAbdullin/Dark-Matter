"""Audit a q64 continuation against the three already archived q32 intervals."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_coupled_serial_restart_reference import main,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_matched_explicit_support_coupled_reference import sources
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def audit(D,eos,work,limit=3):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")
    S=a[0]["snapshot"]
    assert a[0]["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    ix=len(p)-1
    def m(s):
        return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                       [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    n0,_,e0=m(S)
    N0=n0.sum(axis=1)*S["scale_factor"]**3
    E0=e0.sum()*S["scale_factor"]**3
    old=import_restart(S["positive_exterior_restart"]).array()
    output=[]
    for j,v in enumerate(a[1:limit+1],1):
        file=work/"saved"/f"{j}.json"
        if not file.exists() and j==3:
            file=work/"q64"/"checkpoint.json"
        B=read(file)
        assert B["blocks_finished"]==j and B["blocks_total"]==3
        assert B["plan"]["source_sha256"]==sha(work/"source64.json")
        assert B["plan"]["eos_sha256"]==sha(eos)
        R,Z=v["snapshot"],B["snapshot"]
        assert v["snapshot_sha256"]==digest(R) and B["snapshot_sha256"]==digest(Z)
        b=Z["scale_factor"]/S["scale_factor"]
        k=Z["scale_factor"]**3
        assert abs(R["scale_factor"]/Z["scale_factor"]-1)<1e-13
        assert abs(R["T_visible_GeV"]/Z["T_visible_GeV"]-1)<1e-13
        N32,K32,E32=m(R)
        N64,K64,E64=m(Z)
        store=import_restart(Z["positive_exterior_restart"])
        C=store.array()
        check=C[:,:old.shape[1]].copy()
        check[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(check,old,rtol=3e-14,atol=0.)
        H=np.sqrt((TabulatedRadiation(eos,Z["T_visible_GeV"]).state(1.)[1]+E64.sum()+store.energy(1.))/3)/MPL
        hi=np.array([high(Z,1./b,n)[1:3] for n in (8,14)])/(H*E64.sum())
        assert max(abs(hi[1]/hi[0]-1))<.01
        included=C[5,old.shape[1]:].sum()*k
        total=(N64.sum()*k+included)/N0.sum()-1
        assert abs(total)<5e-11
        output.append(dict(block=j,b_relative_to_first_q32=b,T_GeV=Z["T_visible_GeV"],
                           q32_sha256=v["snapshot_sha256"],q64_sha256=B["snapshot_sha256"],
                           shared_count_L1_per_initial=(np.sum(abs(N32-N64[:,:ix]),axis=1)*k/N0).tolist(),
                           shared_full_energy_L1_over_total_initial=(np.sum(abs(E32-E64[:,:ix]),axis=1)*k/E0).tolist(),
                           newly_represented_number_fraction=(N64[:,ix:].sum(axis=1)*k/N0).tolist(),
                           newly_represented_energy_fraction=(E64[:,ix:].sum(axis=1)*k/E0).tolist(),
                           included_combined_number_residual=float(total),
                           high_prompt_inverse_and_tagged_energy_over_current_rho_H=hi.tolist()))
    return dict(scope="Same three expansion intervals and originally archived q32 solution. The "
                      "explicit q64 branch begins with empty padded shells while preserving the "
                      "same older collisionless gravity. All endpoint comparisons conditional on "
                      "the selected coupled model; high H source includes only bare inverse decay.",
                source_q32_sha256=a[0]["snapshot_sha256"],controls=output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("work",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--execute",action="store_true")
    p.add_argument("--limit",type=int,choices=(1,2,3),default=3)
    x=p.parse_args()
    if x.execute:
        sources(x.work,x.tapes)
        a=read(x.tapes/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")
        b=a[-1]["snapshot"]["scale_factor"]/a[0]["snapshot"]["scale_factor"]
        (x.work/"saved").mkdir(exist_ok=True)
        for j in range(1,4):
            backup=x.work/"saved"/f"{j}.json"
            if backup.exists():
                continue
            main(x.eos,x.work/"source64.json",x.work/"q64",b,3,True,1)
            f=x.work/"q64"/"checkpoint.json"
            assert read(f)["blocks_finished"]==j
            backup.write_bytes(f.read_bytes())
    print(json.dumps(audit(x.tapes,x.eos,x.work,x.limit),indent=2,allow_nan=False))
