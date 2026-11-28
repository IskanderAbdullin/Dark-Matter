"""Compare two independent real-absorption charts on the same finite soft band."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_hard_photo_complement_reference import audit
from finite_mass_qed_parent_integrated_signed_soft_reference import Source


def integrate_soft(S,lo,hi,n):
    # No virtual counterterm here; compare only the same real scattering
    # before matching. Source.signed uses unit angular average and k^2
    # scaled squared matrix element relative to the Born spin trace.
    pref=2*((S.M+S.m)**2-S.s**2)/(64*np.pi**3*S.q*S.Es)
    x,w=roots_legendre(n)
    vals=[]
    for u in x:
        k=(lo+hi)/2+(hi-lo)*u/2
        amp=S.signed(-k,42,10,24)[0,:2]
        bose=1/np.expm1(np.longdouble(k))
        vals.append(pref*amp*np.array([bose,1+bose])/k)
    return np.asarray((hi-lo)/2*(w@np.asarray(vals)),dtype=float)


def run():
    fixture=[]
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        lo,hi=4e-6,1e-5
        S=Source(d)
        a,b=(integrate_soft(S,lo,hi,n) for n in (8,14))
        rows=[]
        for mesh in ((7,6,22),(11,10,32)):
            low=audit(d,1.,58.,*mesh,lo)
            high=audit(d,1.,58.,*mesh,hi)
            overlap=[low[k]-high[k] for k in ('gain','loss')]
            rows.append(dict(mesh=mesh,cumulative_lo=[low[k] for k in ('gain','loss')],
                             cumulative_hi=[high[k] for k in ('gain','loss')],finite_band=overlap))
        fixture.append(dict(grid=g,band_edges=[lo,hi],
                            increasing_normed_incoming_and_outgoing_separately=b.tolist(),
                            soft_photon_quadrature_relative_defect=float(max(abs(a/b-1.))),
                            inverse_chart=rows,
                            inverse_chart_refinement=max(abs(np.divide(rows[-1]['finite_band'],rows[-2]['finite_band'])-1.)),
                            chart_relative_defects=(np.asarray(rows[-1]['finite_band'])/b-1.).tolist()))
    return dict(scope='Finite shared band of the real absorption only; no virtual term or new physical rate.',data=fixture)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
