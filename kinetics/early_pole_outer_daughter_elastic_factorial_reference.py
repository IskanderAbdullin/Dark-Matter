"""Attribute the outside daughter effect under matched elastic channel toggles."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def compare(native,outer):
    if set(native)!=set(outer) or set(native)!={'none','unlike','same','both'}:
        raise ValueError('Require all four matched elastic toggles.')
    expected={'none':(0.,0.,'frozen'),'unlike':(1.,0.,'full'),
              'same':(0.,1.,'frozen'),'both':(1.,1.,'full')}
    common=('duration','steps','projection','grid','aperture_grid','hard_edge_sets','compton_strength',
            'conversion_strength','pair','method','response_smear')
    n0, o0=native['none'],outer['none']
    q={}
    functions={}
    for k in expected:
        a,b=native[k],outer[k]
        for n,x in ((n0,a),(o0,b)):
            assert all(n[c]==x[c] for c in common)
            assert (x['strength'],x['same_soft_strength'],x['screening_update'])==expected[k]
            for s in ('hard_moments','cold_bins','instantaneous_source_by_q'):
                np.testing.assert_array_equal(n['snapshots'][0][s],x['snapshots'][0][s])
        for c in common:
            assert a[c]==b[c]
        ca,cb=[conservation(n) for n in (a,b)]
        va,vb=(np.asarray(z['state_history']) for z in (a,b))
        assert va.shape==vb.shape==(a['steps']+1,)
        diff=(vb-vb[0])-(va-va[0])
        functions[k]=diff
        q[k]=dict(native_mean=a['scalar_gain']/a['duration'],outside_mean=b['scalar_gain']/b['duration'],
                  mean_outside_effect=float(diff[-1]/a['duration']),conservation=[ca,cb])
    c=functions['both']-functions['unlike']-functions['same']+functions['none']
    return dict(cases=q,signed_elastic_interaction_in_outside_effect=float(c[-1]/n0['duration']),
                outside_effect_by_macro={k:x.tolist() for k,x in functions.items()},
                scope='Conditional finite count effects of the same outside daughter preparation at separately held channel controls.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',nargs=4,required=True,metavar=('NONE','UNLIKE','SAME','BOTH'))
    p.add_argument('--outer',nargs=4,required=True,metavar=('NONE','UNLIKE','SAME','BOTH'))
    a=p.parse_args()
    def read(paths):
        return {k:json.loads(Path(f).read_text()) for k,f in zip(('none','unlike','same','both'),paths,strict=True)}
    print(json.dumps(compare(read(a.native),read(a.outer)),indent=2,allow_nan=False))
