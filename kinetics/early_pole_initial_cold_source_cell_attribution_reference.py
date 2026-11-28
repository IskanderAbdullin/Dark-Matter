"""Partition the prepared occupation effect by charged hard cell at held dispersion."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel


def run(bundle,order):
    with gzip.open(bundle,'rt') as f:
        d=json.load(f)['branches']
    a,b=(d[p] for p in ('steps','cell-conservative'))
    sa,sb=[x['snapshots'][0] for x in (a,b)]
    N=np.asarray(sa['cold_bins'])
    np.testing.assert_array_equal(N,sb['cold_bins'])
    C=Cold(a['grid'],1,0,len(N),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    M=C.M
    for x,y,w in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
        np.testing.assert_array_equal(x,w)
    E,F,*_=M.weights(np.asarray(sa['hard_moments']),True)
    _,H,*_=M.weights(np.asarray(sb['hard_moments']),True)
    _,_,_,fs,pole=C.scalar(N)
    np.testing.assert_array_equal(F[M.r[0]],H[M.r[0]])

    def calc(z):
        x=np.array([w*band_kernel(C,E,z,q,pole,f,order)
                    for q,w,f in zip(C.q,C.sw,fs,strict=True)])
        np.testing.assert_allclose(x[:,0]-x[:,1],x[:,2],atol=3e-14,rtol=3e-11)
        return x

    base=calc(F)
    np.testing.assert_allclose(base[:,2],sa['instantaneous_source_by_q'],rtol=1e-8,atol=3e-12)
    report={}
    for role,j in (('parent',1),('daughter',2)):
        group=M.r[j]
        e=M.edge_sets[j]
        Z=F.copy()
        Z[group]=H[group]
        direct=calc(Z)-base
        effects=[]
        out=[]
        for k,i in enumerate(range(group.start,group.stop)):
            Y=F.copy()
            Y[i]=H[i]
            dk=calc(Y)-base
            effects.append(dk)
            out.append(dict(hard_momentum_interval=[e[k],e[k+1]],gain_inverse_net=dk[:,:3].sum(axis=0).tolist(),
                            max_logit_offset=float(np.max(abs(H[i]-F[i])))))
        total=np.asarray(effects).sum(axis=0)
        defect=float(np.max(abs(total-direct)))
        assert defect<6e-14
        x=total[:,2].sum()
        t=np.array([a['gain_inverse_net'][2] for a in out])
        q=np.argsort(-abs(t))
        report[role]=dict(signed_direct_gain_inverse_net=direct[:,:3].sum(axis=0).tolist(),
                          additivity_maximum_mode_defect=defect,max_three_cells=[out[i] for i in q[:3]],
                          top_three_absolute_fraction=float(abs(t[q[:3]]).sum()/abs(t).sum()),
                          n_cells_with_nonzero_net_effect=int(np.count_nonzero(abs(t)>1e-15)),
                          n_cells_with_distinct_occupation=int(sum(a['max_logit_offset']>1e-10 for a in out)),
                          signed_sum_relative_to_absolute_of_net_changes=float(x/abs(t).sum()),
                          all_cells=out)
    return dict(order=order,scope='Static conditional probability insertion on fixed native dressed dispersion; no dynamics.',
                aperture=a['aperture_grid'],reports=report)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--order',type=int,default=12)
    a=p.parse_args()
    if a.order<4:
        p.error('Use at least four radial nodes.')
    print(json.dumps(run(a.bundle,a.order),indent=2,allow_nan=False))
