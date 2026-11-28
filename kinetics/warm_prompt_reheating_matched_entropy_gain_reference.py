"""Contrast gross regeneration on passive spectra and matched-count maximum-entropy shapes."""
import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Cells,maxent
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,m,mp,Es
from warm_prompt_reheating_heavy_gain_probe_reference import cells,rates


def run(eos,cache,previous):
    old=json.loads(previous.read_text())
    input_case={(j["TR_GeV"],j["T_GeV"]):j["controls"][-1]["extra_Q_per_H_per_entropy"] for j in old["rows"]}
    b=StandardModel(eos)
    out=[]
    for tr,t in ((.05,.035),(.1,.085),(.1,.075),(.1,.06),(.1,.05),(.1,.035)):
        y=json.loads((cache/f"{tr:.3f}-{t:.3f}-171.json").read_text())
        n=np.r_[y["stable_yields"][1],sum(y["stable_yields"])]
        p=np.asarray(y["scalar_differential_bins"])*t
        rho,_,s,_=b.state(t)
        H=np.sqrt(rho/3)/MPL
        E=np.asarray(y["redshifted_means"])
        target=s*(E[:2]@np.asarray(y["stable_yields"])+E[2]*n[0])
        passive=input_case[tr,t]
        ctl=[]
        for n_grid,ang in ((170,22),(238,32)):
            edge=np.r_[0.,np.geomspace(.0025,90.,n_grid)]*t
            species=[Cells(edge,mp,1),Cells(edge,m,4)]
            T,z=maxent(species,n*s,target,t)
            a=np.array([v[0] for v in z])
            assert abs(a.sum()/target-1)<3e-12
            p0,p1=(cells(edge,k[2]/v/s) for k,v in zip(z,n,strict=True))
            r0,R0=rates(p0[0],p0[1]*s*n[0],mp,ang,True)
            r1,R1=rates(p1[0],p1[1]*s*n[1],m,ang,True)
            total=2*(r0.sum()+r1.sum())/(s*H)
            ctl.append(dict(cells=n_grid,angle_order=ang,dark_T_over_visible=T/t,
                            energy_ledger_relative_error=abs(a.sum()/target-1),
                            scalar_chemical_potential_over_dark_T=z[0][3],
                            stable_chemical_potential_over_dark_T=z[1][3],
                            entropy_reference_extra_Q_per_H_per_entropy=total,
                            rate_over_passive=total/passive,
                            net_scalar_energy_transfer_over_reference_scalar_energy_per_H=float(((Es/M)*R1.sum()-(1-Es/M)*R0.sum())/(H*a[0])),
                            scalar_incoming_fraction=float(r0.sum()/(r0.sum()+r1.sum()))))
        assert abs(ctl[0]["entropy_reference_extra_Q_per_H_per_entropy"]/ctl[1]["entropy_reference_extra_Q_per_H_per_entropy"]-1)<.004
        out.append(dict(TR_GeV=tr,T_GeV=t,input_probe_extra_Q_per_H_per_entropy=passive,controls=ctl))
    return dict(scope="Imposed finite-box BE/FD maximum entropy comparator at the same two numbers and total dark energy. "
                      "Not a proof that this state is reached. Gross Pauli-free regeneration evaluated separately and "
                      "not fed back. No claim about an evolved or thermalized dark sector.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("eos","cache","previous"):
        p.add_argument(v,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.cache,a.previous),indent=2,allow_nan=False))
