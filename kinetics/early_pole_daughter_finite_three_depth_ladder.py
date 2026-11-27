"""Actual local paired gain over nested hard geometries at common cold32 sampling."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


def read(path):
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as z:
            return json.load(z)
    return json.loads(Path(path).read_text())


def one(v, depth, fixture):
    native, local=v
    base=fixture[0]
    use=('projection','duration','report_panels','rtol','base_atol','history_exposure','cold_momentum_nodes_over_T')
    for x,arm in zip(v,('native','local-daughter'),strict=True):
        if (x.get('aperture_subdivision',1)!=depth or x.get('prepared_charged_swap','native')!=arm
                or any(x[k]!=base[k] for k in use) or abs(x['energy_defect'])>=3e-12
                or abs(x['entropy_defect'])>=3e-11 or abs(x['hard_number_defect'])>=3e-12
                or not x.get('elastic_channels',True) or x.get('selected_elastic_channel','both')!='both'):
            raise ValueError(f"Invalid matched scope or ledger at depth {depth}: {arm}.")
    if depth!=1:
        flags=local['local_archive_support']
        if (len(set(flags['indices']))!=2*depth
                or flags['momentum_interval']!=[1.4138888888888888,1.4666666666666666]):
            raise ValueError(f"Nested daughter intervention not reproduced at depth {depth}.")
    t=np.array(native['history_exposure'])
    raw=[]
    for x in v:
        f=np.array(x['cold_mode_count_history'])
        if f.shape!=(len(t),32) or not np.all(np.isfinite(f)&(f>0.)):
            raise ValueError("Unexpected cold history shape or value.")
        np.testing.assert_array_equal(f.T.copy(order='C').sum(axis=0),x['cold_count_history'])
        np.testing.assert_array_equal(f[0],base['cold_mode_count_history'][0])
        raw.append(f)
    mode=(raw[1]-raw[1][0])-(raw[0]-raw[0][0])
    z=np.array([fsum(row) for row in mode])
    return dict(depth=depth,paired_mean=float(z[-1]/t[-1]),
                first_cross_brackets=[t[j:j+2].tolist() for j in range(len(z)-1) if z[j]<0<=z[j+1]],
                terminal_positive=int(np.count_nonzero(mode[-1]>0)),
                terminal_mode_gains=mode[-1].tolist(),paired_counts=z.tolist())


def run(pairs):
    levels=[one(v,d,pairs[0]) for d,v in zip((1,2,4,8)[:len(pairs)],pairs,strict=True)]
    rates=np.array([x['paired_mean'] for x in levels])
    d=np.diff(rates)
    sizes=[int(np.count_nonzero(np.array(w['terminal_mode_gains'])<u['terminal_mode_gains']))
           for u,w in zip(levels[:-1],levels[1:],strict=True)]
    return dict(scope="The actual paired hard ladder at a fixed aperture and time; no continuum bound.",
                exposure=pairs[0][0]['history_exposure'],
                levels=levels,adjacent_mean_shifts=d.tolist(),
                last_shift_over_first=float(d[-1]/d[0]) if d[0] else None,
                modes_smaller_on_each_adjacent_refinement=sizes)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for d in (1,2,4):
        p.add_argument(f'--hard{d}',nargs=2,required=True)
    p.add_argument('--hard8',nargs=2)
    a=p.parse_args()
    pairs=(a.hard1,a.hard2,a.hard4)+((a.hard8,) if a.hard8 is not None else ())
    print(json.dumps(run([[read(k) for k in q] for q in pairs]),indent=2,allow_nan=False))
