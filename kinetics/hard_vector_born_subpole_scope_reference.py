"""Measure additional finite Born support while staying below the vector pole."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from born_thermal_quark_hard_portal_reference import XEDGES
from conditional_relic_normalization_reference import conversion
from disjoint_hard_thermal_measured_source_reference import HardComplement
from measured_partial_hadron_reheat_fold_reference import Fold, HadronFold, spectra
from three_species_decay_onsager_reference import mp
from timelike_portal_onsager_coverage_reference import V


def run(eos, kloe, bes):
    F=Fold(eos,(8,21,16),6,upper=.50)
    H=HadronFold(F,spectra(kloe,bes,10,"linear"),6)
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    def measured(T):
        l,h=H.sample(T)
        q=l[:2]+h[:2,:,0].sum(axis=0)
        return np.array([f*(q[0]+2*q[1]),l[2:].sum()+h[:2,:,1].sum()])
    out=[]
    for tmin,t_probe in ((.25,.3245),(.30,.3368)):
        rr=[]
        for order in (8,12):
            rows=[]
            for x in (14,16,20,24):
                tmax=.345
                assert x*tmax<V
                B=HardComplement(F,tmin,tmax,order,np.r_[XEDGES,np.arange(13.,x+1)])
                def b(T):
                    return B.sample(T)[1]*np.array([f,1.])
                root=brentq(lambda T: measured(T)[0]+b(T)[0]-.1188,.27,.344,xtol=2e-9)
                rows.append(dict(born_x_max=x,max_nominal_mass_GeV=x*tmax,root_GeV=root,
                                 disjoint_counts_and_energy_at_fixed_T=b(t_probe).tolist()))
            assert np.all(np.diff([z["root_GeV"] for z in rows])<0.)
            rr.append(dict(order=order,rows=rows))
        diff=float(max(abs(a["root_GeV"]-b["root_GeV"]) for a,b in zip(rr[0]["rows"],rr[1]["rows"])))
        assert diff<2e-9
        out.append(dict(tmin_GeV=tmin,t_probe_GeV=t_probe,cases=rr,max_order_difference_GeV=diff))
    return dict(controls=out, input_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (eos,kloe,bes)],
           vector_pole_GeV=V,
           scope="Born-only sensitivity for finite below-resonance support on a prescribed clock. Same measured mask; no strict correction in this comparison, complete Born tail, resonance term, physical QCD uncertainty, or comparison to observed abundance.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("eos","kloe","bes"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.bes),indent=2,allow_nan=False))
