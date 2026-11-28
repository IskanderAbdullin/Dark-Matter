"""Factorial channel controls for the local and its outside compensation."""
import argparse
import json
from pathlib import Path

import numpy as np


CHANNELS=('none','unlike','same','both')
SHAPES=('native','local-daughter','local-complement-only')


def compare(paths):
    base=paths['both'][0]
    controls=('projection','duration','report_panels','rtol','base_atol','history_exposure')
    measured={}
    receipts={}
    for ch in CHANNELS:
        x=paths[ch]
        a=paths['both']
        all_counts=[]
        for key,p,q in zip(SHAPES,x,a,strict=True):
            if p.get('prepared_charged_swap','native')!=key or any(p[k]!=base[k] for k in controls):
                raise ValueError('Prepared input or solver controls differ.')
            if p.get('complement_retune')!=q.get('complement_retune'):
                raise ValueError('Complement coefficients changed with channel selection.')
            label=p.get('selected_elastic_channel','both') if p.get('elastic_channels',True) else 'none'
            if label!=ch:
                raise ValueError('Recorded channel does not match requested control.')
            if not (abs(p['energy_defect'])<3e-12 and abs(p['entropy_defect'])<3e-11
                    and abs(p['hard_number_defect'])<3e-12):
                raise ValueError('Invalid terminal energy, entropy or hard count closure.')
            c=np.asarray(p['cold_count_history'])
            if c.shape!=(base['report_panels']+1,) or c[0]!=base['cold_count_history'][0]:
                raise ValueError('Different exposure array or seed.')
            all_counts.append(c-c[0])
        n,l,k=all_counts
        measured[ch]=dict(local=l-n,complement=k-n)
        t=base['duration']
        receipts[ch]=dict(mean_native=float(n[-1]/t),local_mean_effect=float((l[-1]-n[-1])/t),
                     complement_mean_effect=float((k[-1]-n[-1])/t),
                     maximum_energy_defect=max(abs(p['energy_defect']) for p in x),
                     maximum_entropy_defect=max(abs(p['entropy_defect']) for p in x))
    effects={}
    for key in ('local','complement'):
        d={ch:measured[ch][key] for ch in CHANNELS}
        inter=d['both']-d['unlike']-d['same']+d['none']
        effects[key]=dict(final_interaction=float(inter[-1]/base['duration']),
                          sampled_interaction=inter.tolist(),
                          sampled_effect={ch:p.tolist() for ch,p in d.items()})
    return dict(cases=receipts,preparation_channel_interactions=effects,
                exposure=base['history_exposure'],
                scope='Nonadditive response to model channel switches at fixed dielectric and geometry.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in CHANNELS:
        p.add_argument('--'+n,nargs=3,required=True,help='Native, aperture local, complement only.')
    a=p.parse_args()
    print(json.dumps(compare({k:[json.loads(Path(f).read_text()) for f in xs]
                              for k,xs in vars(a).items()}),indent=2,allow_nan=False))
