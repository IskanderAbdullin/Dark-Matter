"""Compare stitched direct FD histories with the instantaneous algebraic lag."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_expansion_forced_modes import parcel
from sm_dark_entropy_exchange_reference import StandardModel


def run(path,start,continuation):
    raw=[x.read_bytes() for x in (start,continuation)]
    d,e=[json.loads(x) for x in raw]
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if any(q['table_input_sha256']!=h for q in (d,e)) or e['archive_sha256']!=hashlib.sha256(raw[0]).hexdigest():
        raise ArithmeticError('Need exact paths, continuation and common EOS')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    out=[]
    for i,(a,b) in enumerate(zip(d['rows'],e['rows'],strict=True)):
        if a['operator']!=b['operator']:
            raise ArithmeticError('Distinct histories cannot be joined')
        u,v=np.array(a['controlled_full_path']),np.array(b['points'])
        if np.max(abs(u[-1]-v[0]))>2e-15:
            raise ArithmeticError('Observed phase-space restart changed bytes')
        rows=[]
        for x,T,eta,_,beta,*_ in np.r_[u,v[1:]]:
            g,k=[parcel(eos,z,n,T)['operators'][i] for z,n in ((p,112),(q,168))]
            t,s=[np.array(j['algebraic_eta_and_beta']) for j in (g,k)]
            if (k['operator']!=a['operator'] or np.max(abs(t-s))>2e-7
                    or np.max(abs(np.array(g['positive_mode_rates_over_H'])/k['positive_mode_rates_over_H']-1.))>2e-6):
                raise ArithmeticError('Local coefficients differ over the same actual trajectory')
            rows.append(dict(x=x,T_SM=T,actual_eta_and_beta=[eta,beta],frozen_algebraic_eta_and_beta=s.tolist(),
                             slow_fast_rate_over_H=k['positive_mode_rates_over_H'],
                             actual_over_algebraic_eta=eta/s[0] if abs(s[0])>1e-12 else None,
                             lag_difference_local_two_mesh_max=np.max(abs(t-s))))
        out.append(dict(operator=a['operator'],rows=rows))
    return dict(rows=out,input_sha256=[hashlib.sha256(x).hexdigest() for x in raw],table_input_sha256=h,
                scope='Only these two finite direct FD histories and local equilibrium linearizations. No time-independent lag, chemical decoupling or final density is inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('eos','start','continuation'):
        p.add_argument(x,type=Path)
    o=p.parse_args()
    print(json.dumps(run(o.eos,o.start,o.continuation),indent=2,allow_nan=False))
