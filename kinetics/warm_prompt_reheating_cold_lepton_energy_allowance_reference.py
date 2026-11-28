"""Add the conservative late charged-lepton energy/count allowance to the cold scalar ceiling."""
import argparse
import json
from pathlib import Path
import numpy as np

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import M,m,mp
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_late_visible_recycling_envelope_reference import COUP2,VECTOR_MASS,energy_rate_coefficient


def run(eos,isolated):
    d=json.loads(isolated.read_text())
    b=Background(eos)
    scale=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    eta=288*COUP2/(np.pi**3*VECTOR_MASS**4)
    ans=[]
    for r in d["anchors"]:
        T=r["visible_T"]
        _,s=visible(b,T)
        H=r["photon_H_floor_at_anchor"]
        A=float(np.exp((b.gs(np.log(T))-np.log(2))/3))
        Tbar=A*T
        C,_=energy_rate_coefficient(Tbar)
        I=C*Tbar**9/(4*H)
        e=eta*Tbar**5/(3*H)
        budget=(r["initial_dark_energy_density"]+I)*np.exp(e)
        N=I/(m*s)
        rec=d["bounded_heavy_coefficient"]*budget**2/(M*M*H*s)
        assert M>=2*m and rec>=r["extra_weighted_yield_upper_if_isolated"]
        ans.append(dict(T_visible_GeV=T,additional_comoving_pair_energy_upper=I,
                        elastic_integrating_factor_exponent=e,remaining_energy_ceiling=budget,
                        direct_weighted_yield_upper=N,regeneration_weighted_yield_upper=rec,
                        total_scalar_omega_upper_including_leptons=scale*(r["initial_weighted_yield"]+N+rec)))
    return dict(scope="Specified late vacuum charged-lepton pair and elastic energy input with a photon-only "
                      "expansion floor and at least two entropy degrees of freedom. Positive inputs are granted "
                      "from the start for the energy ceiling. Does not include late hadron processes, additional "
                      "energy sources or non-vacuum widths; initial anchors remain conditional.",anchors=ans)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path);p.add_argument("isolated",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.isolated),indent=2,allow_nan=False))
