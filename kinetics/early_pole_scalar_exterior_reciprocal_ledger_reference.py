"""Compute matched hard and scalar ledgers for omitted decay source packets."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID,occup


def run(path,projection,order):
    with gzip.open(path,'rt') as f:r=json.load(f)['branches'][projection]
    g=r['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:a=json.load(f)
    q=np.asarray(a['scalar_momentum_edges'])
    fa=np.maximum(0.,np.asarray(a['endpoint_scalar_signed_shell_number_over_T3']))/(np.diff(q**3)/(6*np.pi**2))
    C=Cold(g,1,0,len(r['snapshots'][0]['cold_bins']),True,exposure=r['duration'],
           reference_aperture_grid=r['aperture_grid'],archive_projection=projection)
    for x,y in zip(C.M.edge_sets,r['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    s=r['snapshots'][-1]
    Y=np.asarray(s['hard_moments'])
    E,F,G,*_=C.M.weights(Y,True)
    _,_,_,_,mass=C.scalar(np.asarray(s['cold_bins']))
    assert np.max(abs(E-G))<1e-12
    nodes,weights=roots_legendre(order)
    out=[]
    for lo,hi in ((1e-5,.001),(.001,.1),(.1,5.)):
        ln=np.log(np.unique(np.r_[lo,hi,q[(lo<q)&(q<hi)]]))
        h=np.diff(ln)[:,None]/2
        qs=np.exp(((ln[:-1]+ln[1:])[:,None]/2+h*nodes).ravel())
        ws=(h*weights).ravel()*qs**3/(2*np.pi**2)
        rates=np.zeros(4);hard=np.zeros_like(Y);ent=0.
        for p,w in zip(qs,ws,strict=True):
            fs=float(occup(p,q,fa))
            if fs<=0:
                # Archive emptiness is a real singular entropy edge; do
                # not substitute a floor into the reported entropy.
                raise ValueError(f"Empty scalar occupations are outside this entropy control at {p!r}.")
            v,dy,ds=band_kernel(C,E,F,p,mass,fs,radial_order=12,ledger=True)
            rates+=w*v
            ent+=w*ds
            hard+=w*dy
        N=[float(hard[z,0].sum()) for z in C.M.r[1:]]
        En=[float(np.sum(G[z]*hard[z])) for z in C.M.r[1:]]
        energy=sum(En)+rates[3]
        assert max(abs(N[0]+rates[2]),abs(N[1]-rates[2]))<3e-21
        assert abs(energy)<3e-18 and ent>=0.
        out.append(dict(q_interval=[lo,hi],scalar_net_count_and_energy=rates[-2:].tolist(),
                        parent_daughter_counts=N,parent_daughter_energies=En,
                        number_defect=max(abs(N[0]+rates[2]),abs(N[1]-rates[2])),
                        energy_defect=energy,positive_entropy_rate=ent))
    return dict(grid=g,projection=projection,order=order,snapshot_exposure=s['exposure'],rows=out,
                scope='Frozen reciprocal derivative only; no exterior trajectory or photon scattering invoked.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--order',type=int,default=8)
    a=p.parse_args()
    print(json.dumps(run(a.bundle,a.projection,a.order),indent=2,allow_nan=False))
