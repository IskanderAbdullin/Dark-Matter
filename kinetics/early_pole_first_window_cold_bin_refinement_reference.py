"""Compare joined first-window histories while changing only the scalar cold-bin mesh."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation,select


def compare(runs):
    coarse=runs[0]
    keys=('grid','projection','duration','steps','pair','method','strength','same_soft_strength',
          'compton_strength','conversion_strength','screening_update','response_smear','aperture_grid',
          'hard_edge_sets','elastic_energy_budget')
    nc=len(coarse['snapshots'][0]['cold_bins'])
    zc=np.asarray(coarse['state_history'])
    if nc<2:
        raise ValueError('Use a retained cold mesh with at least two cells.')
    out={}
    for run in runs:
        if any(run[k]!=coarse[k] for k in keys):
            raise ValueError('Only the cold scalar mesh may differ.')
        r0,rt=[run['snapshots'][i] for i in (0,-1)]
        c0,ct=[coarse['snapshots'][i] for i in (0,-1)]
        z=np.asarray(run['state_history'])
        assert z.size==zc.size==run['steps']+1 and abs(z[0]-zc[0])<3e-17
        assert z[-1]-z[0]==run['scalar_gain']
        np.testing.assert_array_equal(r0['hard_moments'],c0['hard_moments'])
        n=len(r0['cold_bins'])
        if n%nc or n<nc or n!=run.get('scalar_bins',8):
            raise ValueError('The finer scalar grid must strictly nest the coarse cells.')
        def row(x,y,s):
            fine=np.asarray(x[s]).reshape(nc,n//nc).sum(axis=1)
            base=np.asarray(y[s])
            d=fine-base
            return dict(coarsened_fine=fine.tolist(),native_coarse=base.tolist(),signed_gap=float(d.sum()),
                        L1_difference=float(np.abs(d).sum()),
                        L1_fraction_of_coarse=float(np.abs(d).sum()/np.abs(base).sum()))
        c_in=row(r0,c0,'cold_bins')
        assert c_in['L1_difference']<4e-17
        s_in=row(r0,c0,'instantaneous_source_by_q')
        s_end=row(rt,ct,'instantaneous_source_by_q')
        gains=(z-z[0])-(zc-zc[0])
        out[str(n)]=dict(energy=conservation(run),initial_counts=c_in,
                         initial_instantaneous_rates=s_in,terminal_instantaneous_rates=s_end,
                         mean_rate=run['scalar_gain']/run['duration'],
                         gain_difference_fraction=(run['scalar_gain']-coarse['scalar_gain'])/coarse['scalar_gain'],
                         largest_absolute_gain_difference=float(np.max(abs(gains))))
    if list(map(int,out))!=sorted(set(map(int,out))):
        raise ValueError('Sort the nested meshes in increasing order.')
    return dict(projection=coarse['projection'],duration=coarse['duration'],steps=coarse['steps'],
                aperture=coarse['aperture_grid'],meshes=out,
                scope='Fixed archived cold support; compare volume-overlap initial occupations and joined evolution only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='+')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(compare([select(x,a.projection) for x in a.paths]),indent=2,allow_nan=False))
