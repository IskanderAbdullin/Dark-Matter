"""Join the verified native fine past bridge without filling the absent future bridge."""
import argparse
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_native_scalar_bridge_join import ledger,field,normalize


def run(low,bridge,cold):
    lo=next(c for c in low if c['grid']=='fine')
    base=next(c for c in cold if c['grid']=='fine')
    if len(bridge)!=1 or bridge[0]['grid']!='fine' or bridge[0]['scope']!=[.0003,.0015] or lo['joined_support']!=[1e-6,.0003]:
        raise ArithmeticError('Need exact adjacent original fine pieces and a separately checked past bridge')
    b=bridge[0]
    raw=np.asarray(base['integrated_number_scalar_energy_and_excess_coefficients'])
    mass=(raw[1,0]-raw[2,0])/raw[0,0]
    omega=base['pole_omega']
    if (np.hypot(.0015,mass)>=omega or b['omega']!=omega or abs(mass/b['m']-1.)>2e-15
        or b['past_or_future_counts']!=[b['panels'],0]
        or not np.array_equal(np.asarray(lo['pieces'][0]['raw']),raw)):
        raise ArithmeticError('The common original frozen pole or the past support changed')
    u=np.sum([normalize(p['raw'],mass) for p in lo['pieces']],axis=0,dtype=np.longdouble)
    v=normalize(b['number_scalar_energy_and_centered_coefficients'],mass)
    q=field(v,omega,mass)[2]
    if np.max(abs(q-b['charged_energy_coefficients']))>2e-13*np.max(abs(q)):
        raise ArithmeticError('Independent charged ledger mismatch')
    return dict(grid='fine',pole_omega=omega,scalar_mass=mass,verified_contiguous_ranges=[[1e-6,.0003],[.0003,.0015]],
                deliberately_open_ranges=[[.0015,.004],[.004,.015]],original_low_only=ledger([u],omega,mass),
                joined_past=ledger([u,v],omega,mass),
                scope='Fine conditional three-sector source on the verified past-support union only; not the future interval or all momentum.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('low','bridge','cold'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    def j(f):
        return json.loads(f.read_text())
    print(json.dumps(run(j(a.low),j(a.bridge),j(a.cold)),indent=2,allow_nan=False))
