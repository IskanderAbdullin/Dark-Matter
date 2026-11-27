"""Conditional pointwise cold Bose gain from an explicit heavy-energy dose."""

import argparse
import json
from pathlib import Path

import numpy as np

from decay_scalar_ir_relaxation_reference import limits
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import Es,M,P,m,mp,width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def envelope(cold,cap,n):
    gamma=cap/M
    assert M<cap<M*Es/mp and 0<cold<P
    boost=np.sqrt((gamma-1)*(gamma+1))
    Es_lb=Es/(gamma+boost)+(mp**2/(Es+P))*boost
    start=np.sqrt((Es_lb-mp)*(Es_lb+mp))
    if start>=cold:
        return dict(cold_reach_floor_GeV=float(start),coefficient_over_rest_width_GeV_minus4=0.,panels=0)
    p=np.geomspace(start,cold,n+1)
    E,L,_=limits_vec(p)
    # In this p-range E_min decreases as p grows. Endpoint inequalities
    # give a true upper over each band, not just a sample of the curve.
    assert np.all(np.diff(L)<0.) and L[-1]>M
    np.testing.assert_allclose(L[0],cap,rtol=8e-15)
    pref=lambda x,e:4*M*M/(2*P*x*e)
    heavy=lambda l:(2/np.pi**2)*l*l*np.sqrt((l-M)*(l+M))
    upper=pref(p[:-1],E[:-1])/heavy(L[1:])
    at=pref(p,E)/heavy(L)
    assert np.max(upper)>=np.max(at)
    return dict(cold_reach_floor_GeV=float(start),coefficient_over_rest_width_GeV_minus4=float(np.max(upper)),
                sampled_coefficient_GeV_minus4=float(np.max(at)),panels=n)


def limits_vec(p):
    return np.array([limits(float(q)) for q in p]).T


def run(eos,source,reactivated):
    x=read(source)
    S=x["snapshot"]
    z=read(reactivated)
    assert digest(S)==x["snapshot_sha256"]==z["base_sha256"]
    T=S["T_visible_GeV"]
    _,_,E=get_mom(np.asarray(S["lab_momentum_edges_GeV"]),
                   [S["entropy_variable_"+k] for k in ("scalar","fermion")])
    rho=E.sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
    r=[]
    for pad in (1.05,2.,3.5,7.5,20.):
        Ecap=float(np.hypot(pad,mp)+np.hypot(pad,m))
        out=[envelope(.1*T,Ecap,n) for n in (24,96,384)]
        if out[-1]["panels"]:
            assert out[0]["coefficient_over_rest_width_GeV_minus4"]>=out[1]["coefficient_over_rest_width_GeV_minus4"]>=out[2]["coefficient_over_rest_width_GeV_minus4"]
        d=float(out[-1]["coefficient_over_rest_width_GeV_minus4"])*width*rho/H
        t=[dict(hard=r["hard_momentum_GeV"],
                fixed_benchmark_for_parent_energy_fraction=r["eligible_stationary_heavy_energy_over_rho"],
                hypothetic_uniform_exponent_over_short_step_using_benchmark_as_cap=float(d*1e-7*r["eligible_stationary_heavy_energy_over_rho"]))
           for r in z["rows"] if r["eligible_stationary_heavy_energy_over_rho"]>0.] if pad==7.5 else []
        r.append(dict(common_daughter_cap_GeV=pad,implied_parent_energy_cap_GeV=Ecap,refinement=out,
                max_exponent_per_H_time_if_hard_energy_equals_anchor_dark_energy=d,
                required_integrated_hard_energy_fraction_times_Hdt_for_one_percent_Bose_factor=(float(np.log1p(.01)/d) if d else None),
                illustrative_not_verified_through_time_benchmarks=t))
    return dict(scope="Upper comparison for prompt gain of ln(1+f_s) with a specified finite parent cap and "
             "integrated eligible-heavy energy. Other terms, redshift and actual time-dependent energy "
             "are excluded. Frozen energy benchmarks are not verified as time-dependent bounds.",cases=r)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("source",type=Path)
    p.add_argument("reactivated",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.source,a.reactivated),indent=2,allow_nan=False))
