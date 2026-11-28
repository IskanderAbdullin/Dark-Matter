"""Retain one common hard support while evolving a manufactured parent/daughter input swap."""
import argparse
import json
from unittest.mock import patch

import numpy as np

import early_pole_cold_screened_split_reference as split
from early_pole_prepared_hard_swap_reference import read,retune


def run(bundle,case,steps,no_elastic=False,only_unlike=False,only_same=False,hold_screening=False):
    if sum((no_elastic,only_unlike,only_same))>1:
        raise ValueError("Elastic controls are mutually exclusive.")
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
            if case.startswith(('local-','outside-')):
                cut=(1.4138888888888888,1.4666666666666666)
                role=(1,2) if case.endswith('-both') else ((1,) if case.endswith('-parent') else (2,))
                outside=case.startswith('outside-')
                boxes={}
                for j in role:
                    edge=np.asarray(M.edge_sets[j])
                    sel=np.flatnonzero((edge[:-1]>=cut[0]-2e-14)&(edge[1:]<=cut[1]+2e-14))
                    if len(sel)!=2 or not np.allclose(edge[[sel[0],sel[1]+1]],cut,rtol=0.,atol=2e-14):
                        raise ValueError("The chosen two hard cells are not present on this aperture.")
                    if outside:
                        sel=np.setdiff1d(np.arange(len(edge)-1),sel,assume_unique=True)
                    ix=M.r[j].start+sel
                    y[ix]=fine_cont[ix]
                    boxes[str(j)]=np.c_[edge[sel],edge[sel+1]].tolist()
                metadata['outside_aperture_swaps' if outside else 'localized_aperture_swaps']=boxes
            # Build all subsequent channels and their reference response
            # relative to this exact prepared state, not a pre-swap state.
            M.init=y
    with patch.object(split,'Screened',PreparedInput):
        r=split.run(a['grid'],a['duration']/a['steps'],steps,True,aperture_grid=a['aperture_grid'],
                    strength=0. if no_elastic or only_same else 1.,
                    screening_update='frozen' if no_elastic or only_same or hold_screening else 'full',
                    same_soft_strength=0. if no_elastic or only_unlike else 1.,compton_strength=1.,conversion_strength=1.,
                    snapshots=True,elastic_energy_budget=None if no_elastic else 5e-15)
    r.update(prepared_charged_swap=case,manufactured_hard_shape_control=True,**metadata)
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--case',choices=('native','parent','daughter','both','matched-parent','matched-daughter','matched-both','local-parent','local-daughter','local-both','outside-parent','outside-daughter','outside-both'),required=True)
    p.add_argument('--steps',type=int,default=16)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--no-elastic',action='store_true')
    group.add_argument('--only-unlike',action='store_true')
    group.add_argument('--only-same',action='store_true')
    p.add_argument('--hold-screening',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
