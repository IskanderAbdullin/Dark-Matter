"""Frozen cold-current stress with moment-preserving high exterior reactivation."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M, mp, width
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_late_cold_parent_time_dilation_bound_reference import vacuum_coefficients
from warm_late_cold_scalar_prompt_partition_reference import angular_counts, recoil_gate
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_late_cold_scalar_single_use_receipt_reference import run as receipt


def run(eos,anchor):
    x=read(anchor)
    S=x["snapshot"]
    assert x["snapshot_sha256"]==digest(S)
    T=S["T_visible_GeV"]
    cutoff=.1*T
    cut=np.hypot(mp,cutoff)
    Eh_min=recoil_gate(np.array([cutoff]),1.)[0][0]
    ph_min=np.sqrt((Eh_min-M)*(Eh_min+M))
    reference=receipt([anchor])["cases"][0]["controls"][1]
    coef=vacuum_coefficients(cutoff)[1]["vacuum_rate_coefficient_upper"]
    base_f=reference["archived_cold_BE_max_f"]
    results=[]
    for hard in (.713,1.55,3.3,6.5):
        D,Er,err,*_=reactivate(S,hard,.02,cap=7.5)
        p=np.asarray(D["lab_momentum_edges_GeV"])
        N,_,E=get_mom(p,[D["entropy_variable_"+k] for k in ("scalar","fermion")])
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/3.)/MPL
        norm=N[0].sum()*H
        r=[]
        for closed in (True,False):
            check=[angular_counts(D,np.array([cut]),order,closed)/norm for order in (6,10,14)]
            ref=check[-1][0,:2].max()
            if ref>0:
                assert max(abs(check[-1][0]-check[-2][0]))<2e-4*ref
            else:
                assert all(np.all(a[0]==0.) for a in check)
            r.append(dict(closed=closed,orders=[6,10,14],rates=[a.tolist() for a in check]))
        # The reactivation does not touch the retained cold cells.
        oldp=np.asarray(S["lab_momentum_edges_GeV"])
        k=int(np.searchsorted(oldp,cutoff))
        assert np.array_equal(p[:k+1],oldp[:k+1])
        assert np.array_equal(np.asarray(D["entropy_variable_scalar"])[:,:k],
                              np.asarray(S["entropy_variable_scalar"])[:,:k])
        V=np.array([high(D,float(ph_min),order) for order in (8,12)])
        allowance=(V[-1,3]/M)*(width/norm)*(1+base_f)*coef
        gain=r[-1]["rates"][-1][0][0]
        assert gain<=allowance*(1+2e-10)
        if V[-1,3]>0:
            assert abs(V[-2,3]/V[-1,3]-1.)<1e-4
        results.append(dict(eligible_stationary_heavy_energy_over_rho=float(V[-1,3]/E.sum()),
                            eligible_gross_inverse_count_over_Ns_H=float(V[-1,0]/norm),
                            cold_gross_decay_count_cap_over_Ns_H=float(allowance),
                            cold_gross_decay_count_fraction_of_bound=(float(gain/allowance) if allowance>0 else None),
                            hard_momentum_GeV=hard,halfwidth=.02,shape_sha256=digest(D),
                            moment_reconstruction_error=err,reactivated_energy_over_rho=float(Er/E.sum()),
                            cold_gate_p_over_T=.1,controls=r))
    return dict(scope="Frozen prompt current including intentionally injected hard shapes made from preserved receipts. "
                      "No new source, claimed observed spectrum, repeat accounting or cold-population trajectory.",
                base_sha256=x["snapshot_sha256"],rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("anchor",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.anchor),indent=2,allow_nan=False))
