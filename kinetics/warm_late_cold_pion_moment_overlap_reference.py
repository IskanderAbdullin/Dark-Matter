"""Conditionally subtract the measured overlap before spending two positive moments."""
import argparse,hashlib,json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre
from kloe_pion_band_portal_reference import table,NB_PER_GEV2
from warm_babar_dispersive_gap_reference import mu_kernel_lower,MU
from warm_light_occupation_stage_audit_reference import read


def kloe_moments(grid,cs,err,last,multiples,n=16):
    # Conditional piecewise-constant spectral envelope in each KLOE bin.
    hi=np.minimum(grid+.005,last)
    lo=grid-.005
    take=hi>lo
    hi,lo=hi[take],lo[take]
    w=np.maximum(cs[take]-multiples*err[take],0.)/NB_PER_GEV2
    edges=lo[:,None]+(hi-lo)[:,None]*np.linspace(0,1,n+1)
    mu=mu_kernel_lower(edges[:,1:].ravel()).reshape(len(lo),n)/(4*np.pi**3)
    upper_s_analytic=float(w@np.sum(np.diff(edges,axis=1)*mu,axis=1))
    # Check the lower against an independent smooth integral, not against
    # statistical evidence for its piecewise-constant cross section.
    x,wx=roots_legendre(32)
    y=(x+1)/2
    q=(lo[:,None]+(hi-lo)[:,None]*y).ravel()
    H=((y*y*(1-y))*wx/2)@(1/(y*y+(1-y)*q[:,None]/MU**2)).T
    numer=float(w@((H.reshape(len(lo),len(y))@(wx/2))*(hi-lo))/(4*np.pi**3))
    assert 0<upper_s_analytic<=numer
    return upper_s_analytic,numer


def run(gap,kloe):
    x=read(gap)
    gd,cs,err=table(kloe)
    U,V=x["bands"]
    assert U["s_GeV2"][1]==V["s_GeV2"][0] and V["s_GeV2"][1]==.1
    a,b=[v["controls"][-1]["max_yield_per_dispersion_unit"] for v in (U,V)]
    prior_C,prior_B=[v["assumed_positive_dispersion_cap"] for v in (U,V)]
    assert b>a
    out=[]
    for s in (0,1,2,3):
        a5,chk5=kloe_moments(gd,cs,err,.50**2,s)
        a6,chk6=kloe_moments(gd,cs,err,.63**2,s)
        C=prior_C-a6; B=prior_B-a5
        if C<0. or B<0.:
            out.append(dict(marginal_subtractions=s,compatible=False,remaining_prior_ceiling=[C,B]))
            continue
        # Coefficients are nonnegative, so for the tightest total upper
        # first spend as much shared moment as allowed in the bigger ratio.
        high=min(C,B); low=C-high
        upper=a*low+b*high
        assert upper<=U["controls"][-1]["cold_yield_upper"]+V["controls"][-1]["cold_yield_upper"]
        out.append(dict(marginal_subtractions=s,compatible=True,
            measured_moment_lower_to_point50=a5,measured_moment_lower_to_point63=a6,
            independent_exact_integral_to_point50=chk5,independent_exact_integral_to_point63=chk6,
            shared_remaining_dispersion_cap=C,separate_gap_dispersion_cap=B,
            extremal_budget_split_threshold_gap=[low,high],cold_yield_joint_upper=upper,
            held_fraction_upper=upper/x["held_equality_cold_yield"]))
    return dict(scope="Illustrative deterministic consistency of whole theory and BABAR upper inputs with "
                "piecewise KLOE lower spectral assumptions on their overlapping measured ranges. "
                "No statistical combination, covariance inference, or claim that the assumed lower spectrum is covered.",
                KLOE_sha256=hashlib.sha256(kloe.read_bytes()).hexdigest(),scenarios=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("gap",type=Path);p.add_argument("kloe",type=Path)
    a=p.parse_args();print(json.dumps(run(a.gap,a.kloe),indent=2,allow_nan=False))
