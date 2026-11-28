"""Integrate frozen-snapshot tail diagnostics without evolving the tail."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator


def integrate_samples(raw):
    data=raw['cases']
    t=np.array([x['exposure'] for x in data])
    if len(t)<4 or not np.all(np.diff(t)>0.) or t[0]!=0.:
        raise ValueError('Use the ordered --snapshots all output.')
    support=[x['interval'] for x in data[0]['rows']]
    if any([r['interval'] for r in x['rows']]!=support for x in data):
        raise ValueError('Snapshot band bounds do not match.')
    y=np.array([[(r['net'],r['scalar_energy_source']) for r in x['rows']] for x in data])
    x=trapezoid(y,x=t,axis=0)
    z=PchipInterpolator(t,y,axis=0).integrate(t[0],t[-1])
    ii=support.index([1e-6,1e-5])
    rows=[]
    for j,m in enumerate(support):
        rows.append(dict(interval=m,trapezoid_net_count_and_energy=x[j].tolist(),
                         pchip_net_count_and_energy=z[j].tolist(),
                         pchip_minus_trapezoid=(z[j]-x[j]).tolist()))
    return dict(grid=raw['grid'],projection=raw['projection'],sample_times=t.tolist(),panels=rows,
                exterior_relative_to_cold=dict(trap=((x.sum(0)-x[ii])/x[ii]).tolist(),
                                               pchip=((z.sum(0)-z[ii])/z[ii]).tolist()),
                scope='Interpolation of frozen source snapshots only, not backreacting abundance or a bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input')
    a=p.parse_args()
    print(json.dumps(integrate_samples(json.loads(Path(a.input).read_text())),indent=2,allow_nan=False))
