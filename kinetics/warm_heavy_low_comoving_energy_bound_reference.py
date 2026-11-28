"""Conditional heavy energy envelope in finite comoving bands."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import M,mp,P,Es


def read(path):
    return json.loads(Path(path).read_bytes())


def upper(E,delta,source,cap,n,scalar_share=1.):
    """Initially empty heavy modes; bounds below a cut only, not the total."""
    assert cap < M*P/mp  # E_phi,min decreases until scalar can be at rest.
    edge=np.linspace(0.,cap,n+1)
    lower,upper=edge[:-1],edge[1:]
    EH=np.hypot(M,upper)
    # Stable rational energy of the anticollinear scalar. The lower
    # endpoint is decreasing in the entire declared heavy momentum range.
    ehi=(EH*Es+upper*P)/M
    elo=(Es**2+(mp*upper/M)**2)/ehi
    assert np.all(np.diff(elo)<0.)
    klo=np.sqrt((elo-mp)*(elo+mp))
    assert np.min(klo)>0.
    with np.errstate(divide="ignore"):
        # rho_phi b^3<=scalar_share*E. The lab angle introduces b/P_c,
        # so the inferred inverse numerator can only shrink as b^-2.
        A=delta*np.pi**2*M*(scalar_share*E)/(P*lower*klo*elo**2)
    f=np.minimum(delta,(A+source)/(1-delta))
    return float((2/np.pi**2)*np.sum((upper**3-lower**3)*EH*f/3)/E)


def run(pions,barrier,light):
    X,Y,Z=map(read,(pions,barrier,light))
    rows=[]
    for a,b,c in zip(X["cases"],Y["results"],Z["cases"],strict=True):
        assert a["q"]==b["q"]==c["q"]
        E=a["controls"][-1]["new_available_energy_upper_GeV4"]
        h=c["diagnostic_heavy_energy_fraction_for_direct_dose_half"]
        cases=[]
        for k in b["cases"]:
            delta=k["cap"]
            S=delta*(1-delta)*k["total_bounded_input_over_outgoing_decay_margin"]
            points=[]
            for cut in (.001,.003,.01,.03,.10,.30,.60,1.,1.4,2.):
                vals=[upper(E,delta,S,cut,n) for n in (128,512,2048)]
                assert vals[0]>vals[-1]
                points.append(dict(comoving_cut_GeV=cut,
                                   below_cut_energy_fraction_upper=vals[-1],
                                   below_cut_fraction_upper_clipped_by_total=min(1.,vals[-1]),
                                   can_bound_below_cut_to_under_light_gain_budget=bool(vals[-1]<h),
                                   provision_for_above_cut_if_only_direct_light_gains=max(0.,h-vals[-1])))
            cases.append(dict(light_cap=delta,checked_subbands=points))
        rows.append(dict(q=a["q"],direct_light_diagnostic_heavy_fraction=h,cases=cases))
    return dict(scope="Heavy modes initially empty; no heavy elastic redistribution. "
                      "Assumes a common uniform light cap, earlier E budget and all-mode "
                      "number-changing input bounds. Finite scalar energy but no scalar "
                      "occupancy cap. The complementary high-comoving-momentum band "
                      "is not depleted or bounded by this low-band estimate.",
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                              for p in (pions,barrier,light)],cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pions","barrier","light"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.pions,a.barrier,a.light),indent=2,allow_nan=False))
