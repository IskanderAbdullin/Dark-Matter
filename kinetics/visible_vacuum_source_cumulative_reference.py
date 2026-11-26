"""Independent cumulative count in the prescribed, vacuum source limit."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from empty_lepton_pair_spectrum_reference import Source
from three_species_decay_onsager_reference import M,m
from sm_dark_entropy_exchange_reference import StandardModel, MPL


def accumulate(eos, TR, TF, order=(9,24,16), time_order=8):
    """Integrate Y = n/s directly in visible ln(T), with no FRW momentum grid."""
    s_light, s_heavy = [Source(order,dark_mass=k) for k in (m,M)]
    x,w=roots_legendre(time_order)
    lo,hi=np.log([TF,TR])
    edge=np.unique(np.r_[lo,np.clip(eos.gs.x,lo,hi),hi])
    result=np.zeros((2,2))
    sinit=eos.state(TR)[2]
    for a,b in zip(edge[:-1],edge[1:]):
        for l,weight in zip((a+b)/2+(b-a)/2*x,(b-a)/2*w):
            T=np.exp(l)
            rho,_,s,cv=eos.state(T)
            H=np.sqrt(rho/3)/MPL
            # c_V = T ds/dT; d ln a = -c_V/(3s) d ln T.
            factor=weight*cv/(3*H*s)
            a=(sinit/s)**(1/3)
            for i,src in enumerate((s_light,s_heavy)):
                q,_,_,n,_=src.event_moments(T)
                result[i]+=factor*np.array([2*n.sum()/s,a**4*np.sum(q*n)])
    return result


def test(eos, raw=None):
    TR=.05
    ends=np.array([.045,.04,.035,.03,.025,.02,.015,.012])
    full=np.array([accumulate(eos,TR,t) for t in ends])
    n,energy=full[...,0],full[...,1]
    fine=accumulate(eos,TR,.012,order=(12,32,20),time_order=12)
    assert np.all(np.diff(n,axis=0)>=0.)
    max_order=float(np.max(abs(full[-1]/fine-1.)))
    assert max_order<3e-3
    val=None
    if raw is not None:
        previous=np.array(raw["stable_yields"])
        val=float(np.max(abs(n[0]/previous-1.)))
        assert val<3e-5
    return dict(reheating_T_GeV=TR,
                baseline_spectrum_relative_crosscheck=val,
                high_order_cumulative_relative_difference=max_order,
                channels=["direct light stable", "heavy prompt stable or scalar (equal production)"],
                rows=[dict(T_visible_GeV=float(t),yield_per_channel=y.tolist(),
                           net_injected_a4_energy_per_pair_channel=e.tolist(),
                           fraction_of_12_MeV_cumulative=(y/n[-1]).tolist(),
                           fraction_of_12_MeV_cumulative_injected_energy=(e/energy[-1]).tolist())
                      for t,y,e in zip(ends,n,energy)])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("--spectra")
    a=p.parse_args()
    r=None
    if a.spectra:
        r=next(x for x in json.loads(Path(a.spectra).read_text())["rows"]
               if x["TR"]==.05 and abs(x["TF"]-.045)<1e-8)
    print(json.dumps(test(StandardModel(a.eos),r),indent=2,allow_nan=False))
