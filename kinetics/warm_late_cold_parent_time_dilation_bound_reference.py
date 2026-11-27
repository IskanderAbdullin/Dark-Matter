"""Bound the open cold-decay rate of preexisting hard parents by count or energy."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import Es, M, P, m, mp, width
from warm_late_cold_scalar_pauli_count_bound_reference import max_filled_length

from warm_late_cold_scalar_single_use_receipt_reference import run as recorded
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def vacuum_coefficients(q):
    C=np.hypot(mp,q)
    ymin=np.arcsinh(P/mp)-np.arcsinh(q/mp)
    ymax=np.arcsinh(P/mp)+np.arcsinh(q/mp)
    bmin=np.tanh(ymin)
    lo,hi=np.cosh(ymin),np.cosh(ymax)
    y=np.linspace(ymin,ymax,8001)
    gamma=np.cosh(y)
    angular=np.maximum(C-mp*np.cosh(y-np.arcsinh(P/mp)),0.)/(2*P*np.sinh(y))
    v=[]
    for r in (1,2):
        # E_{s,min} >= Es/(2 gamma), sqrt(gamma^2-1)>=bmin*gamma.
        # The analytic comparison curve may be maximized on a larger
        # domain than the actual opening, which can only raise its maximum.
        n=r+1
        where=(n+1)*Es/(2*n*C)
        assert lo<where<hi
        upper=C/((n+1)*2*P*bmin*where**n)
        sampled=float(np.max(angular/gamma**r))
        assert sampled<=upper*(1+3e-14) and sampled>0.98*upper
        v.append(dict(power_of_parent_gamma=r,vacuum_rate_coefficient_upper=float(upper),
                      sampled_vacuum_max=float(sampled),analytic_upper_gamma=float(where)))
    return v


def run(eos,paths):
    out=[]
    for raw,old in zip(paths,recorded(paths)["cases"],strict=True):
        x=read(raw)
        S=x["snapshot"]
        assert old["hash"]==x["snapshot_sha256"]==digest(S)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        T=S["T_visible_GeV"]
        N,_,E=get_mom(p,[S["entropy_variable_"+s] for s in ("scalar","fermion")])
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/3.)/MPL
        ext=import_restart(S["positive_exterior_restart"]).array()
        light=N[1].sum()+ext[5,ext[1]==m].sum()
        cases=[]
        for o in old["controls"]:
            q=o["q_over_T"]*T
            E_min=o["E_parent_gate_GeV"]
            boosted=vacuum_coefficients(q)
            W_min=2*P*np.sqrt((E_min/M-1)*(E_min/M+1))
            # Rate = Gamma/gamma times angular measure averaged over
            # its *vacuum* interval, not the branching fraction within
            # the interval divided by an unknown blocked normalization.
            stat=1+o["archived_cold_BE_max_f"]
            for r in boosted:
                r["Bose_adjusted_rate_coefficient_upper"]=stat*r["vacuum_rate_coefficient_upper"]
                r["coefficient_times_rest_width_over_H"]=r["Bose_adjusted_rate_coefficient_upper"]*(width/H)
            cases.append(dict(p_over_T=o["q_over_T"],parent_gate_GeV=E_min,
                              vacuum_interval_width_min_GeV=float(W_min),
                              optional_known_count_blocking_measure_over_min_width=float(max_filled_length(light)/W_min),
                              rate_bounds=boosted))
        out.append(dict(archive=x["snapshot_sha256"],visible_T_GeV=T,H_GeV=float(H),
                        parent_rest_width_over_H=float(width/H), controls=cases))
    return dict(scope="Frozen scalar occupations in cold band, full open unpolarized decay and any valid light FD shape. "
            "Rate bounds multiply Gamma by eligible parent number (r=1) or eligible parent energy / M (r=2). "
            "No bound on the heavy population, production history or changing Bose occupations.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archives),indent=2,allow_nan=False))
