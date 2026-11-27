"""Retain one common hard support while evolving a manufactured parent/daughter input swap."""
import argparse
import json
from unittest.mock import patch

import numpy as np

import early_pole_cold_screened_split_reference as split
from early_pole_prepared_hard_swap_reference import read,retune


def run(bundle,case,steps):
    raw=read(bundle)
    a,b=[raw[p] for p in ('steps','cell-conservative')]
    sa,sb=[q['snapshots'][0] for q in (a,b)]
    fine_steps=np.asarray(sa['hard_moments'])
    fine_cont=np.asarray(sb['hard_moments'])
    for x,y in zip(a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    np.testing.assert_array_equal(sa['cold_bins'],sb['cold_bins'])
    metadata={}
    original=split.Screened
    class PreparedInput(original):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            M=self.M
            np.testing.assert_array_equal(M.init,fine_steps)
            np.testing.assert_array_equal(self.ns,sa['cold_bins'])
            for x,y in zip(M.edge_sets,a['hard_edge_sets'],strict=True):
                np.testing.assert_array_equal(x,y)
            p,d=M.r[1:]
            y=M.init.copy()
            if case in ('parent','matched-parent','both','matched-both'):
                y[p]=fine_cont[p]
            if case in ('matched-parent','matched-both'):
                y,check=retune(M,y,M.init,p)
                metadata['parent_moment_match']=check
            if case in ('daughter','both','matched-daughter','matched-both'):
                y[d]=fine_cont[d]
            if case in ('matched-daughter','matched-both'):
                y,check=retune(M,y,M.init,d)
                metadata['daughter_moment_match']=check
            # Build all subsequent channels and their reference response
            # relative to this exact prepared state, not a pre-swap state.
            M.init=y
    with patch.object(split,'Screened',PreparedInput):
        r=split.run(a['grid'],a['duration']/a['steps'],steps,True,aperture_grid=a['aperture_grid'],
                    screening_update='full',same_soft_strength=1.,compton_strength=1.,conversion_strength=1.,
                    snapshots=True,elastic_energy_budget=5e-15)
    r.update(prepared_charged_swap=case,manufactured_hard_shape_control=True,**metadata)
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--case',choices=('native','parent','daughter','both','matched-parent','matched-daughter','matched-both'),required=True)
    p.add_argument('--steps',type=int,default=16)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
