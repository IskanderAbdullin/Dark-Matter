"""Check the early thermal heat timescale against continuing FD pair production."""
import argparse
import json
from pathlib import Path
import numpy as np

from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source,spectrum
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,m,mp
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_signed_slow_closure_reference import closed,quantum_state
from warm_prompt_reheating_quantum_heat_gap_reference import thermal,capacity


def run(eos):
    eos0=StandardModel(eos)
    b=Background(eos)
    D=[]
    for offset in (.5,1.,1.5,2.,2.5,3.,5.,10.,20.,50.,100.):
        tv=.1-offset*1e-6
        r=spectrum(eos0,.1,tv,order=(12,35,24),nout=14)
        rv,sv=visible(b,tv)
        Y=np.r_[sum(r["stable_yields"]),r["stable_yields"][1]]
        L=np.array(r["redshifted_means"])
        E=L[:2]@r["stable_yields"]+L[2]*Y[1]
        T,_,z=closed(Y*sv,sv*E,tv)
        T,eta,N,_,_=quantum_state(Y*sv,sv*E,T,z,14,True)
        H=np.sqrt((rv+sv*E)/3)/MPL
        l=Source((12,35,25),tail=72.,dark_mass=m).event_moments(tv)[3].sum()
        h=Source((12,35,25),tail=72.,dark_mass=M).event_moments(tv)[3].sum()
        influx=2*np.array([l+h,h])/(H*sv)
        K=thermal(T,eta,86,54)
        c=np.array([capacity(T,eta[0],m,4,False,24),capacity(T,eta[1],mp,1,True,24)])
        gap=K*np.sum(1/c)/H
        clock=float(max(influx/Y))
        D.append(dict(drop_below_reheat_keV=offset,dark_T_GeV=T,initial_yields=Y.tolist(),
                      heat_gap_over_H=gap,max_fractional_source_over_H=clock,
                      heat_gap_to_source=gap/clock))
    rates=np.asarray([a["heat_gap_to_source"] for a in D])
    assert np.all(np.diff(rates)>0.) and rates[0]<1.<rates[-1]
    i=np.searchsorted(rates,1.)
    return dict(scope="Early diagnostic on independent inherited direct FD source moments matched to two fixed-number thermal species. "
                      "Linear resonant heat gap versus fractional source growth, not proof of kinetic equilibration or a "
                      "freeze-in transition. Source slopes normalized to the imposed passive baseline; no dark feedback before each sample.",
                unity_ratio_bracket_keV=[D[i-1]["drop_below_reheat_keV"],D[i]["drop_below_reheat_keV"]],rows=D)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos),indent=2,allow_nan=False))
