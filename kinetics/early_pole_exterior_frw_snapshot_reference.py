"""Check the reciprocal warm-band kernel on frozen comoving hard states at imposed scale."""
import argparse
import gzip
import json
from types import MethodType

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_joined_screened_frw_reference import exterior
from early_pole_live_joint_hard_packets_frw_reference import Expanding
from early_pole_native_aperture_witness_reference import GRID


def run(bundle,projection):
    with gzip.open(bundle,'rt') as f:
        raw=json.load(f)['branches'][projection]
    grid=raw['grid']
    shot=raw['snapshots'][-1]
    Y=np.array(shot['hard_moments'])
    N=np.array(shot['cold_bins'])
    C=Cold(grid,1,0,len(N),True,exposure=raw['duration'],reference_aperture_grid=raw['aperture_grid'],
           archive_projection=projection)
    M=C.M
    for a,b in zip(M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
    M.weights=MethodType(Expanding.weights,M)
    M.spectators=exterior(M,grid)
    M.eta=0.
    assert max(abs(M.spectators(0.)[0]-M.background))<3e-12
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        archived=json.load(f)
    bounds=np.asarray(archived['scalar_momentum_edges'])
    keep=(C.spq>=.1)&(C.spq<=5.)&(C.spect>0.)
    Q=C.spq[keep]
    V=np.diff(bounds**3)[keep]/(6*np.pi**2)
    Fw=C.spect[keep]/V
    masses,pole=C.scalar(N)[0],C.scalar(N)[-1]
    sample=[]
    for eta in (0.,2e-5,7.2e-5,3e-4,1e-3):
        M.eta=eta
        E,F,G,*_=M.weights(Y,True)
        assert np.max(abs(E-G))<1e-12
        rhs=C.rhs(np.r_[Y.ravel(),N,0.,0.],False)
        orig=rhs[C.j0:C.j0+C.n]
        cr=np.zeros(C.n)
        H0=np.zeros_like(Y)
        for i,(q,v,f) in enumerate(zip(C.q,C.sw,N/C.sw,strict=True)):
            qv,hard,_=band_kernel(C,E,F,q,pole,f,12,ledger=True,scale_log=eta)
            cr[i]=v*qv[2]
            H0+=v*hard
        assert np.max(abs(cr-orig))<1e-10
        assert abs(np.sum(G*H0)+masses@cr)<5e-16
        warm=np.zeros(4)
        hard=np.zeros_like(Y)
        entropy=0.
        unscaled=0.
        for q,v,f in zip(Q,V,Fw,strict=True):
            r,d,s=band_kernel(C,E,F,q,pole,f,12,ledger=True,scale_log=eta)
            warm+=v*r
            hard+=v*d
            entropy+=v*s
            unscaled+=v*band_kernel(C,E,F,q,pole,f,12)[2]
        de=float(np.sum(G*hard)+warm[-1])
        dc=[float(hard[a,0].sum()) for a in M.r[1:]]
        assert abs(de)<2e-20
        assert max(abs(dc[0]+warm[2]),abs(dc[1]-warm[2]))<2e-20 and entropy>0.
        sample.append(dict(log_scale=eta,cold_replay_max_mode_defect=float(np.max(abs(cr-orig))),
                           cold_total_rate=float(cr.sum()),warm_rates=warm.tolist(),
                           hard_count_reciprocity=dc,energy_exchange_defect=de,
                           unscaled_warm_rate_on_same_hard_energy=unscaled,
                           warm_entropy_rate=entropy))
    return dict(projection=projection,grid=grid,history_exposure=raw['duration'],rows=sample,
                scope='Prescribed snapshot scale only; held comoving distributions and no live FRW history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
