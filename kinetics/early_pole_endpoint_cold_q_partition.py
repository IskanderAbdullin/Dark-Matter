"""Compare normalized q profiles and locally signed flux on exactly shared cold edges."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_receipt import read,reduce
from early_pole_endpoint_peer_charge_shape_transport import quantile_support_bounds


def geometry(path):
    x=read(path)
    E=np.array([x[0]['bounds'][0]]+[r['bounds'][1] for r in x])
    z=np.array([r['coeff'] for r in x])
    return E,z


def cumulative(e,w):
    p=np.array(w)/fsum(w)
    return e,np.r_[0.,np.cumsum(p)]


def regroup(e,w,edges):
    ids=np.searchsorted(e,edges)
    if len(edges)<2 or np.any(ids>=len(e)) or not np.array_equal(e[ids],edges) or np.any(np.diff(ids)<=0):
        raise ArithmeticError('Common edges must be exact and ordered in each archive')
    out=np.array([fsum(w[i:j]) for i,j in zip(ids[:-1],ids[1:],strict=True)])
    if abs(fsum(out)/fsum(w)-1)>2e-15:
        raise ArithmeticError('Coarsening lost the number or energy coefficient')
    return out


def shape(E,F,e,f,edges):
    if min(np.min(F),np.min(f))<0 or min(fsum(F),fsum(f))<=0:
        raise ArithmeticError('Need positive finite panel weights')
    a=regroup(E,F,edges)/fsum(F)
    b=regroup(e,f,edges)/fsum(f)
    total=fsum(abs(b-a))/2
    c=np.cumsum(b-a)
    j=np.argsort(-abs(b-a))[:5]
    lo,hi,_=quantile_support_bounds(cumulative(E,F),cumulative(e,f))
    return dict(TV_lower_bound_from_exact_common_bins=total,
                max_shared_boundary_CDF_difference=float(max(abs(c))),
                Wasserstein1_q_panel_supported_interval=[lo,hi],
                important_bins=[dict(q=edges[i:i+2].tolist(),old_new_normalized_mass=[a[i],b[i]],
                                    signed_normalized_difference=b[i]-a[i]) for i in j])


def neutral(e,z,edges):
    B=z[:,0,1]
    D=z[:,0,2]
    n=fsum(B)/fsum(D)
    gb,gd=[regroup(e,v,edges) for v in (B,D)]
    net=(np.longdouble(n)*gd.astype(np.longdouble)-gb.astype(np.longdouble))
    net=np.asarray(net,dtype=float)/fsum(B)
    residual=fsum(net)
    if abs(residual)>5e-15 or abs(fsum(net[net>0])+fsum(net[net<0]))>5e-15:
        raise ArithmeticError('The q-resolved signed redistribution does not sum to zero')
    j=int(np.argmax(abs(net)))
    return dict(number_neutral_occupation=n,positive_flux_sum_over_total_gross_inverse=fsum(net[net>0]),
                negative_flux_sum_over_total_gross_inverse=fsum(net[net<0]),
                net_flux_normalized_by_total_gross_inverse=net.tolist(),
                sum_roundoff=residual,strongest_bin_q=edges[j:j+2].tolist())


def single(a,b):
    (E,Z),(e,z)=map(geometry,(a,b))
    edges=np.intersect1d(E,e)
    if len(edges)<5 or [edges[0],edges[-1]]!=[1e-6,1.5e-5]:
        raise ArithmeticError('Insufficient exact shared boundaries')
    return dict(exact_common_q_edges=edges.tolist(),
        gross_forward_number=shape(E,Z[:,0,0],e,z[:,0,0],edges),
        gross_inverse_number=shape(E,Z[:,0,1],e,z[:,0,1],edges),
        gross_forward_excess=shape(E,Z[:,2,0],e,z[:,2,0],edges),
        neutral_flux_by_chart=[neutral(E,Z,edges),neutral(e,z,edges)])


def run(a,b,levels):
    if tuple(levels) not in ((1,2),(2,4),(4,8),(8,16)):
        raise ArithmeticError('Match consecutive archived scalar levels only')
    data=[]
    for k,files in zip(levels,(a,b),strict=True):
        (a0,A,t),(u,U,s)=map(reduce,files)
        if A!=U or t.keys()!=s.keys() or any(t[j]!=s[j] for j in t if j!='mesh') or t.get('scalar_div',1)!=k or t['grid']!='fine':
            raise ArithmeticError('Panel pair and scalar charts do not agree')
        delta=np.max(abs(a0/u-1))
        if delta>6e-5 or t['time_x']!=.05:
            raise ArithmeticError('Use independently accepted full cold panel integrals')
        data.append(t)
    if data[0].get('charge_div',2)!=2 or data[1].get('charge_div')!=2:
        raise ArithmeticError('Changing the charged chart is outside this test')
    low,high=[single(a[j],b[j]) for j in (0,1)]
    if low['exact_common_q_edges']!=high['exact_common_q_edges']:
        raise ArithmeticError('Coarsened q partition changed under integration')
    compare=('gross_forward_number','gross_inverse_number','gross_forward_excess')
    at='TV_lower_bound_from_exact_common_bins'
    err=[abs(low[c][at]-high[c][at]) for c in compare]
    keys=('positive_flux_sum_over_total_gross_inverse','negative_flux_sum_over_total_gross_inverse')
    flux=max(abs(r[k]-s[k]) for r,s in zip(low['neutral_flux_by_chart'],high['neutral_flux_by_chart'],strict=True) for k in keys)
    if max(*err,flux)>6e-5:
        raise ArithmeticError('The selected shape and subband net controls did not close between meshes')
    return dict(levels=levels,pole_scope='same numerical charge partition, independently fitted pole per scalar chart',
                **high,absolute_pair_mesh_difference_TV=err,absolute_pair_mesh_difference_neutral_redistribution=flux,
                scope='Finite positive panel measures; coarsened TV and W1 enclosures do not interpolate within panels or claim real momentum/distribution limits.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left',type=Path,nargs=2,required=True)
    p.add_argument('--right',type=Path,nargs=2,required=True)
    p.add_argument('--levels',type=int,nargs=2,required=True)
    z=p.parse_args()
    print(json.dumps(run(z.left,z.right,z.levels),indent=2,allow_nan=False))
