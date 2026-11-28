"""Prescribed-background size of the omitted vacuum production after 12 MeV."""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import rule
from empty_lepton_pair_spectrum_reference import Source
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,M


def interval(bg,lo,hi,order,time):
    src=[Source(order,dark_mass=x) for x in (m,M)]
    x,w=rule(time)
    a,b=np.log([lo,hi])
    edges=np.unique(np.r_[a,np.clip(bg.sg.x,a,b),b])
    count=np.zeros(2)
    for l,r in zip(edges,edges[1:]):
        for u,weight in zip((l+r)/2+(r-l)/2*x,w*(r-l)/2):
            T=np.exp(u)
            s=2*np.pi**2/45*np.exp(bg.entropy_log(u))
            H=np.sqrt(np.pi**2/90*np.exp(4*u+bg.rg(u)))/MPL
            measure=weight*(3+bg.sg.derivative()(u))/(3*H*s)
            for j,f in enumerate(src):
                count[j]+=measure*2*f.event_moments(T)[3].sum()
    return count


def test(file,base):
    bg=TabulatedRadiation(file,.012)
    rows=[]
    den=np.array(base["yield_per_channel"])
    for lo,hi in ((.0093,.012),(.004,.0093),(.001,.004)):
        a=interval(bg,lo,hi,(8,18,14),5)
        b=interval(bg,lo,hi,(11,28,18),8)
        d=np.max(abs(a/b-1.))
        assert d<3e-4
        rows.append(dict(T_GeV=[lo,hi],additional_yield_per_channel=b.tolist(),
                         relative_to_each_12_MeV_cumulative_yield=(b/den).tolist(),
                         relative_quadrature_change=float(d)))
    return rows


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("radiation_table")
    p.add_argument("source_reference")
    a=p.parse_args()
    raw=json.loads(Path(a.source_reference).read_text())["rows"]
    r=next(x for x in raw if x["T_visible_GeV"]==.012)
    print(json.dumps(dict(rows=test(a.radiation_table,r),scope="prescribed adiabatic visible radiation and vacuum production, not a late inverse-reaction bound"),indent=2,allow_nan=False))
