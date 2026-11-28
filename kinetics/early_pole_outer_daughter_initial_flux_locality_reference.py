"""Project elastic initial directions onto cold-active and outside hard moment cells."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        old=json.load(f)['branches']
    a,b=(old[k] for k in ('steps','cell-conservative'))
    y0,y1=(np.asarray(d['snapshots'][0]['hard_moments']) for d in (a,b))
    n=np.asarray(a['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(n,b['snapshots'][0]['cold_bins'])
    C=Screened(a['grid'],1,len(n),a['duration'],(1,1,2),.08,1.,a['aperture_grid'],'steps')
    M=C.M
    np.testing.assert_array_equal(y0,M.init)
    for x,y,z in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y); np.testing.assert_array_equal(x,z)
    elastic=ElasticChannels(C,a['grid'],1.)
    active=np.zeros(len(y0),dtype=bool)
    daughter=np.zeros_like(active)
    daughter[M.r[2]]=True
    for j,span in enumerate(M.r[1:],1):
        edges=np.asarray(M.edge_sets[j])
        sel=np.flatnonzero((edges[:-1]>=1.4138888888888888-2e-14)&(edges[1:]<=1.4666666666666666+2e-14))
        assert len(sel)==2
        active[span.start+sel]=True
    def source(y):
        return Cold.rhs(C,np.r_[y.ravel(),n,0.,0.],False)[C.j0:-2]
    out={}
    for case,selector in [('native',np.zeros_like(active)),('outer_daughter',daughter&~active)]:
        Y=np.where(selector[:,None],y1,y0)
        _,F,G,*_=M.weights(Y,True)
        U=C.collision(G,F,G)[0]
        V=elastic.soft.collision(G,F,G)[0]
        directions={}
        for k,D in [('unlike',U),('same',V),('both',U+V)]:
            points=[]
            for dt in (3e-15,1e-15):
                d=lambda z:(source(Y+dt*z)-source(Y-dt*z))/(2*dt)
                tot,local,far=[d(x) for x in (D,D*active[:,None],D*(~active)[:,None])]
                size=float(np.max(abs(tot)))
                defect=float(np.max(abs(tot-local-far))/size)
                points.append(dict(probe=dt,total=float(tot.sum()),in_aperture=float(local.sum()),
                                   outside=float(far.sum()),mode_additivity_relative_defect=defect))
            directions[k]=points
        out[case]=dict(initial_rate=float(source(Y).sum()),directions=directions)
    return dict(scope='Initial held-graph derivatives for native and outer-daughter hard inputs. Recompute dressed energies at each trial, not a time integral.',
                channels=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
