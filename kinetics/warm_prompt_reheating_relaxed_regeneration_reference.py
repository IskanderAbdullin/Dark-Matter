"""Re-evaluate regeneration before and after isolated prompt-return relaxation."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Cells
from sm_dark_entropy_exchange_reference import MPL,StandardModel
from three_species_decay_onsager_reference import m,mp
from warm_prompt_reheating_heavy_gain_probe_reference import rates


def run(file,eos):
    d=json.loads(file.read_text())
    r=d["terminal_state"]
    edge=np.asarray(r["momentum_edges_GeV"])
    cells=(Cells(edge,mp,1),Cells(edge,m,4))
    start=np.asarray(r["initial_species_number_and_kinetic_energy"])
    finish=np.asarray(r["species_number_and_kinetic_energy"])
    assert np.max(abs(start[:,0].sum(axis=1)/finish[:,0].sum(axis=1)-1))<3e-12
    shape0=[c.reconstruct(v)[0] for c,v in zip(cells,start,strict=True)]
    shape1=np.asarray(r["species_entropy_variables"])
    b=StandardModel(eos)
    rho,_,entropy,_=b.state(d["T"])
    H=np.sqrt(rho/3)/MPL
    x=[]
    for radial,angle in ((5,17),(8,27)):
        xx,w=roots_legendre(radial)
        mid,half=(edge[:-1]+edge[1:])/2,np.diff(edge)/2
        p=mid[:,None]+half[:,None]*xx
        s=[]
        for iso,want in ((shape0,start),(shape1,finish)):
            cases=[]
            for c,v,q in zip(cells,iso,want,strict=True):
                idx=np.arange(len(c.vol))[:,None]
                E=np.hypot(p,c.mass)
                N=c.g*w*half[:,None]*p*p*c.f(v,E,idx)/(2*np.pi**2)
                assert abs(N.sum()/q[0].sum()-1.)<2e-5
                keep=N.ravel()>0
                cases.append(rates(p.ravel()[keep],N.ravel()[keep],c.mass,angle))
            s.append((2*np.r_[cases[0],cases[1]]/(H*entropy)).tolist())
        x.append(dict(radial=radial,incoming_angle=angle,scalar_opposite_same_extra_Q_per_entropy_H=s,
                      total_ratio_after_to_before=sum(s[1])/sum(s[0])))
    return dict(scope="Gross regeneration on before/after shapes only; no ongoing regeneration in the fast prompt-return integration.",
                short_H_time=d["log_time_H_units"],resolved_relaxation=d["fractional_progress_towards_maximum_entropy"],tiers=x)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("snapshot",type=Path)
    p.add_argument("eos",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.snapshot,a.eos),indent=2,allow_nan=False))
