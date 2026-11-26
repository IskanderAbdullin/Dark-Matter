"""Measure the Bose enhancement on actual archived warm outgoing pair geometry."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import gross,snapshot_functions
from three_species_decay_onsager_reference import m,mp


def snapshot(path):
    raw=Path(path).read_bytes()
    z=json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)
    return z["snapshot"] if "snapshot" in z else z["snapshots"][-1]


def scalar_shell_extrema(s):
    slope,tilt = np.asarray(s["entropy_variable_scalar"])
    assert np.all(slope>0) and np.all(slope+tilt>0)
    f0,f1=1/np.expm1(slope),1/np.expm1(slope+tilt)
    return np.r_[0.,np.minimum(f0,f1),0.],np.r_[0.,np.maximum(f0,f1),0.]


def one(s,order):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    f,_,g,_=snapshot_functions(s)
    _,vac,projection=gross(p,m,mp,g,np.zeros_like,False,*order,project=True)
    _,bose,_=gross(p,m,mp,g,f,False,*order)
    lo,hi=scalar_shell_extrema(s)
    counts=projection[1,0]
    assert len(lo)==len(counts)
    maximum=max(hi)
    lower=(counts@lo)/vac[0]
    upper=(counts@hi)/vac[0]*(1+.5*maximum)
    E=bose[1]/vac[0]-1
    assert lower-2e-15<=E<=upper+2e-15
    pct_exterior=(counts[0]+counts[-1])/(2*vac[0])
    return dict(order=order,vacuum_event_rate_GeV4=float(vac[0]),
                retained_peak=float(maximum),
                retained_number_fraction_in_outgoing_daughters=float(1-pct_exterior),
                exact_positive_Bose_enhancement=float(E),
                shell_enclosure=[float(lower),float(upper)],
                uniform_enhancement_ceiling=float(maximum*(2+maximum)),
                energy_weighted_relative_enhancement=float(bose[2]/vac[2]-1.))


def run(paths):
    data=[]
    for path in paths:
        s=snapshot(path)
        checks=[one(s,order) for order in ((3,7,8,3),(4,10,12,4),(5,13,16,5))]
        a,b=checks[-2:]
        E1,E2=[r["exact_positive_Bose_enhancement"] for r in (a,b)]
        rel=abs(E1/E2-1.) if E2!=0. else (0. if E1==0. else None)
        data.append(dict(path=Path(path).name,T_GeV=s["T_visible_GeV"],
                         sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                         quadratures=checks,last_two_bose_enhancement_relative_difference=rel))
    return dict(scope="Frozen retained warm momentum shapes, no extrapolated exterior "
                      "occupations, forward chi2-antichi2-to-scalar pair geometry only. "
                      "Numerical quadratures; shell enclosures apply within each "
                      "discretely integrated geometry. Not a uniform-in-time Bose cap "
                      "or a bound for a dynamically redistributed reservoir.",
                cases=data)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths",nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.paths),indent=2,allow_nan=False))
