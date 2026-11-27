"""Held per-species momentum kicks from below-pole elastic scattering."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_two_temperature_contact_probe_reference import fold,kgeom


def run(eos,path):
    bg=Background(eos)
    obj=json.loads(path.read_text())
    r=next(p for p in obj["controls"][-1]["heat_arms"] if p["prompt_heat"])["rows"]
    ans=[]
    save=[]
    for i in (4,8,16):
        u=r[i];t=u["visible_T"];Tl,Ts=u["light_and_scalar_dark_T"]
        v,s=visible(bg,t);H=np.sqrt((v+sum(u["energy_per_visible_entropy"])*s)/3)/MPL
        pk=u["light_and_scalar_phase_peak"]
        ps=[]
        for nn in (8,12):
            d=[]
            for j,(T,ma,g) in enumerate(((Tl,m,4),(Ts,mp,1))):
                p,E,w=kgeom(T,ma,nn)
                z=pk[j]*np.exp(-(E-ma)/T)
                w*=g*z/(1+z if j==0 else 1-z)
                assert abs(w.sum()/s/u["yields"][j]-1)<2e-6
                d.append(w@(p*p))
            Q,N,out,K=fold(Tl,Ts,pk,(nn,nn+3,nn+5),exact=True,transport=True)
            ps.append(dict(radial=nn,gross_events_over_H_and_count=(N/(H*s*np.array(u["yields"]))).tolist(),
                           mean_kick_squared_over_H_population_second_moment=(K/(H*np.array(d))).tolist()))
        ans.append(dict(T_visible=t,controls=ps))
    return dict(scope="Positive squared momentum kicks in the single below-pole channel. This moment diagnostic "
                      "does not evaluate an angular relaxation spectrum or prove internal Fermi/Bose closure.",rows=ans)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("handoff",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.handoff),indent=2,allow_nan=False))
