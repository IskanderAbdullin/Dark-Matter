"""Measure the represented prompt high band at exact corrected checkpoints."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M, m, mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def check(eos, prior, paths):
    old=read(prior)["branches"][-1]
    t0=old["T_GeV"]
    b0=old["entropy_expansion_from_original_anchor"]
    T_ref=TabulatedRadiation(eos,t0)
    out=[]
    for path in paths:
        u=read(path)
        S=u["snapshot"]
        assert u["snapshot_sha256"]==digest(S)
        T=S["T_visible_GeV"]
        R=TabulatedRadiation(eos,T)
        a=b0*np.exp((T_ref.s0-R.s0)/3)
        cut=1./a
        p=np.array(S["lab_momentum_edges_GeV"])
        A=[S["entropy_variable_"+k] for k in ("scalar","fermion")]
        E=get_mom(p,A)[2].sum()
        old_density=import_restart(S["positive_exterior_restart"]).energy(1.)
        H=np.sqrt((R.state(1.)[1]+E+old_density)/3)/MPL
        rate=np.array([high(S,cut,n) for n in (6,10,14)])
        rate[:,:3]*=np.array([M,1.,1.])/(E*H)
        rate[:,3]/=E
        support_max=np.hypot(m,p[-1])+np.hypot(mp,p[-1])
        assert np.all(rate>=0) and np.isfinite(rate).all()
        if support_max<=np.hypot(M,cut):
            assert np.all(rate==0.)
            mode="Geometrically closed on represented support."
        else:
            assert np.all(rate>0.)
            mode="Open on represented support."
        out.append(dict(snapshot_sha256=u["snapshot_sha256"],endpoint_T_GeV=T,
                        momentum_cut_GeV=cut,represented_max_pair_energy_GeV=support_max,
                        closure=mode,orders=[6,10,14],normalized_rates=rate.tolist(),
                        fine_relative_rate_shift=None if rate[-1,1]==0. else (rate[-1]/rate[-2]-1).tolist()))
    return dict(scope="Instantaneous resolved support of the selected archived dark phase only; "
                      "quasi-static occupied heavy parent used as the blocking convention. "
                      "No old exterior interactions, coherent heavy evolution or unknown high tails.",
                columns=["high_inverse_count_times_rest_mass_over_rhoH","high_inverse_energy_over_rhoH",
                         "direct_tagged_daughter_energy_over_rhoH","high_quasistatic_heavy_energy_over_rho"],
                cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("prior",type=Path)
    p.add_argument("archive",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(check(a.eos,a.prior,a.archive),indent=2,allow_nan=False))
