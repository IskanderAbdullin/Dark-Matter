"""Compare elastic-centered interleaving at unchanged total elastic work and width."""
import argparse
import json
from unittest.mock import patch

import early_pole_cold_screened_split_reference as split
from early_pole_local_daughter_elastic_subcycle_reference import run as subcycle


def run(bundle,case,macro_steps=16):
    if macro_steps not in (16,32):
        raise ValueError('Require the fixed 64-solve inner comparison.')
    original=split.run
    def alternate(*args,**kwargs):
        return original(*args,**kwargs,outer_order='elastic-slow-elastic')
    with patch.object(split,'run',alternate):
        z=subcycle(bundle,case,macro_steps,32//macro_steps)
    z['elastic_substeps_per_macro_step']=64//macro_steps
    assert sum(len(x['elastic_substeps']) for x in z['elastic'])==64
    return z


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--case',required=True,choices=('native','local-daughter'))
    p.add_argument('--macro-steps',type=int,default=16)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
