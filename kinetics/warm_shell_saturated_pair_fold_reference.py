"""Frozen scalar pair rate on matched shell moments with different sub-shell shapes."""

import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import gross, snapshot_functions, m, mp
from warm_light_occupation_stage_audit_reference import digest, read
from warm_shell_occupation_feasibility_reference import independent_moments, saturated_witness, volume

ORDERS=((2,5,5,2),(3,8,8,3),(4,10,10,4),(6,14,12,5))


def fold(S, p, j, ne):
    a,b=p[j:j+2]
    v=volume(a,b)
    w=saturated_witness(a,b,m,*ne[j])
    lo,hi=np.cbrt(a**3+3*v*(w["inner_volume_fraction"]+np.array([0.,w["interval_volume_fraction"]])))
    assert a<lo<hi<b
    refined=np.sort(np.r_[p,lo,hi])
    Ea,Eb,El,Eh=np.hypot([a,b,lo,hi],m)
    B,_,F,_=snapshot_functions(S)

    def saturated(e):
        return np.where((e>=Ea)&(e<Eb),((e>=El)&(e<Eh)).astype(float),F(e))

    rows=[]
    for r,t,z,k in ORDERS:
        a=[]
        for f in (F,saturated):
            n,h,_=gross(refined,m,mp,f,B,False,r,t,z,k)
            _,g,_=gross(refined,mp,m,B,f,True,r,t,z,k)
            assert h[1]>g[1]>0
            a.append(dict(light_number=n,forward_event_rate=h[1],inverse_event_rate=g[1],
                          net_event_rate=h[1]-g[1],forward_energy_rate=h[2],inverse_energy_rate=g[2]))
        ratios={key:a[1][key]/a[0][key] for key in a[0]}
        rows.append(dict(order=[r,t,z,k],archive=a[0],saturated=a[1],saturated_over_archive=ratios))
    assert abs(rows[-1]["saturated_over_archive"]["light_number"]-1.)<1e-9
    return dict(shell_index=j,shell_number_fraction=ne[j,0]/ne[:,0].sum(),
                interval_GeV=[lo,hi],saturated_witness=w,controls=rows)


def run(directory):
    raw=read(directory/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    S=raw["snapshot"]
    assert raw["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    theta=np.asarray(S["entropy_variable_fermion"])
    ne=np.array([independent_moments(a,b,m,t) for a,b,t in zip(p[:-1],p[1:],theta.T,strict=True)])
    feasibility=read(directory/"warm-shell-occupation-feasibility-2026-11-26.json")
    case=feasibility["cases"][1]["results"][0]
    assert case["snapshot_sha256"]==raw["snapshot_sha256"]
    first=fold(S,p,case["shell_index"],ne)
    second=fold(S,p,int(np.argmax(ne[:,0])),ne)
    assert first["shell_index"]!=second["shell_index"]
    return dict(T_visible_GeV=S["T_visible_GeV"],snapshot_sha256=raw["snapshot_sha256"],
                scope="Frozen retained q32 light-antilight to scalars and inverse only. "
                      "Two separate per-shell replacements with identical radial number and "
                      "kinetic energy; neither is a predicted physical state.",
                cap_audit_shell=first,largest_number_shell=second)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(parser.parse_args().tapes),indent=2,allow_nan=False))
