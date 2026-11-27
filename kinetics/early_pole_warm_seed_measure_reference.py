"""Separate spontaneous quadrature and occupation effects at the same hard state and pole."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID


def run(path, projection, outer_orders, radial_orders):
    with gzip.open(path,'rt') as fh:
        raw=json.load(fh)['branches'][projection]
    g=raw['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as fh:
        spectrum=json.load(fh)
    shot=raw['snapshots'][-1]
    n=np.array(shot['cold_bins'])
    C=Cold(g,1,0,len(n),True,exposure=raw['duration'],
           reference_aperture_grid=raw['aperture_grid'],archive_projection=projection)
    for x,y in zip(C.M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    E,F,*_=C.M.weights(np.array(shot['hard_moments']),True)
    pole=C.scalar(n)[-1]
    take=(C.spq>=.1)&(C.spq<=5.)&(C.spect>0.)
    markers=C.spq[take]
    counts=C.spect[take]
    lo,hi=np.array(spectrum['scalar_momentum_edges'])[:-1][take]**3,np.array(
        spectrum['scalar_momentum_edges'])[1:][take]**3
    vol=(hi-lo)/(6*np.pi**2)
    rows=[]
    for radial in radial_orders:
        matrices={}
        def rates(q):
            if q not in matrices:
                s,_,_,_=band_kernel(C,E,F,q,pole,0.,radial)
                g,i,net,_=band_kernel(C,E,F,q,pole,1.,radial)
                assert abs(g-2*s)+abs(g-i-net)<5e-18
                matrices[q]=np.array([s,i])
            return matrices[q]
        def integrate(q,w,occ,label,order):
            S,I=np.array([rates(x) for x in q]).T
            en=np.hypot(q,pole)
            spont=w*S
            stimulated=spont*occ
            inv=w*I*occ
            net=spont+stimulated-inv
            return dict(label=label,outer_order=order,maximum_occupation=float(occ.max()),
                        spontaneous=float(spont.sum()),stimulated=float(stimulated.sum()),
                        inverse=float(inv.sum()),net=float(net.sum()),
                        spontaneous_scalar_energy=float(en@spont),net_scalar_energy=float(en@net))
        a=integrate(markers,vol,counts/vol,'native_marker',1)
        out=[a]
        decompositions=[]
        for order in outer_orders:
            z,w=roots_legendre(order)
            q=np.cbrt((lo+hi)[:,None]/2+(hi-lo)[:,None]*z/2)
            weights=vol[:,None]*w/2
            E0=np.hypot(markers[:,None],pole)
            En=np.hypot(q,pole)
            de=(q*q-markers[:,None]**2)/(En+E0)
            scaled=de/np.ptp(En,axis=1)[:,None]
            prob=[]
            for y in scaled:
                def p(beta):
                    v=beta*y
                    r=w*np.exp(v-v.max())
                    return r/r.sum()
                beta=brentq(lambda v:p(v)@y,-200.,200.)
                prob.append(p(beta))
            pop=counts[:,None]*np.array(prob)
            assert max(abs((pop*de).sum(axis=1)))<1e-21
            u=integrate(q.ravel(),weights.ravel(),np.broadcast_to((counts/vol)[:,None],q.shape).ravel(),
                        'resolved_uniform_occupation',order)
            b=integrate(q.ravel(),weights.ravel(),(pop/weights).ravel(),'resolved_energy_matched_occupation',order)
            assert u['spontaneous']==b['spontaneous']
            out.extend((u,b))
            dn=a['net']-b['net']
            ds=a['spontaneous']-b['spontaneous']
            do=(a['stimulated']-a['inverse'])-(b['stimulated']-b['inverse'])
            assert abs(dn-ds-do)<1e-21
            decompositions.append(dict(outer_order=order,marker_minus_energy_matched_net=dn,
                spontaneous_difference=ds,occupation_difference=do,spontaneous_fraction_of_net_gap=ds/dn,
                uniform_to_energy_matched_on_identical_nodes=b['net']-u['net']))
        rows.append(dict(radial_order=radial,rows=out,marker_gap_decomposition=decompositions))
    return dict(grid=g,projection=projection,hard_snapshot_exposure=shot['exposure'],fixed_scalar_pole=pole,
                native_marker_nodes=markers.tolist(),shell_bounds=np.c_[np.cbrt(lo),np.cbrt(hi)].tolist(),
                comparisons=rows,scope='Matched fixed-state warm decay rates; no live trajectory or continuum claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--outer-orders',type=int,nargs='+',default=[4,8,12])
    p.add_argument('--radial-orders',type=int,nargs='+',default=[8,16])
    a=p.parse_args()
    if min(a.outer_orders)<3 or min(a.radial_orders)<2:
        p.error('Outer order must be at least three and radial order at least two.')
    print(json.dumps(run(a.bundle,a.projection,a.outer_orders,a.radial_orders),indent=2,allow_nan=False))
