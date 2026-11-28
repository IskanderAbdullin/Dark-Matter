"""Compare two non-nested whole cold bands with native-panel CDF enclosures."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_continuum_cold_band_pilot import run as topology
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def enclosure(edges,mass,cuts):
    edges,mass,cuts=map(lambda v:np.asarray(v,dtype=float),(edges,mass,cuts))
    if (len(edges)!=len(mass)+1 or np.any(np.diff(edges)<=0)
            or not np.all(np.isfinite(mass)) or min(mass)<0
            or not fsum(mass)>0. or min(cuts)<edges[0] or max(cuts)>edges[-1]):
        raise ArithmeticError('Need ordered native edges and a finite positive band')
    cumulative=np.r_[0.,np.cumsum(mass/fsum(mass))]
    completed=np.searchsorted(edges,cuts,side='right')-1
    possibly=np.searchsorted(edges,cuts,side='left')
    if abs(cumulative[-1]-1.)>3e-14:
        raise ArithmeticError('CDF normalization failed')
    return np.clip(cumulative[np.c_[completed,possibly]],0.,1.)


def overlap_tv(e,f,a,b):
    """Bracket the minimum unbinned TV permitted by the two native panel ledgers."""
    q=np.union1d(e,f)
    i,j=[np.searchsorted(z,q[:-1],side='right')-1 for z in (e,f)]
    cap=np.r_[np.asarray(a,dtype=float)/fsum(a),np.asarray(b,dtype=float)/fsum(b)]
    n,m=len(a),len(b)
    C=np.zeros((n+m,len(q)-1))
    ix=np.arange(C.shape[1])
    C[i,ix]=1.
    C[n+j,ix]=1.
    # Common mass can be assigned only inside intervals belonging to one
    # panel in each chart. Two capacities constrain each shared allocation.
    z=linprog(-np.ones(C.shape[1]),A_ub=C,b_ub=cap,bounds=(0,None),
              method='highs',options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    if not z.success:
        raise ArithmeticError('Native overlap solver failed')
    g=np.maximum(z.x,0.)
    for row,c in zip(C.astype(bool),cap,strict=True):
        used=fsum(g[row])
        if used>c:
            g[row]*=c/used
    prices=np.maximum(-z.ineqlin.marginals,0.)
    guard=64*np.finfo(float).eps
    deficit=max(0.,float(max(1.-C.T@prices)))
    prices[:n]+=deficit+guard
    if max(C@g-cap)>guard or min(C.T@prices)<1.:
        raise ArithmeticError('Independent overlap/capacity certificates do not close')
    best=max(0.,1.-fsum(prices*cap)-guard)
    witness=min(1.,1.-fsum(g)+guard)
    if witness<best or witness-best>2e-9:
        raise ArithmeticError('Native interval overlap has no sufficiently tight certificate')
    return dict(minimum_TV_bounds=[best,witness],common_measure_by_union_bin=g.tolist(),
                native_panel_capacity_prices=[prices[:n].tolist(),prices[n:].tolist()])


def compare(e,f,a,b):
    if [e[0],e[-1]]!=[f[0],f[-1]]:
        raise ArithmeticError('Need the same exact retained band')
    q=np.union1d(e,f)
    A,B=enclosure(e,a,q),enclosure(f,b,q)
    separation=np.maximum(0.,np.maximum(A[:,0]-B[:,1],B[:,0]-A[:,1])-64*np.finfo(float).eps)
    i=int(np.argmax(separation))
    overlap=overlap_tv(e,f,a,b)
    if overlap['minimum_TV_bounds'][1]+2e-14<separation[i]:
        raise ArithmeticError('Finite overlap witness contradicts cumulative lower bound')
    return dict(upper_over_lower_band_response=fsum(b)/fsum(a),native_union_cut_q=q.tolist(),
                cdf_enclosures_at_cuts=[A.tolist(),B.tolist()],
                guaranteed_cdf_gap_by_cut=separation.tolist(),
                normalized_TV_lower_bound_via_cdf=float(separation[i]),
                strongest_q=float(q[i]) if separation[i]>0 else None,sharp_compatible_TV=overlap)


def audit(pair,state,archive,grid):
    card,sha=read(state)
    if (card['grid'],card['x_end'],card['subcells_per_archived_cell'],card['steps'])!=(grid,.05,2,2048):
        raise ArithmeticError('Expected native archived D2 endpoint on the same time mesh')
    check=receipt(pair,archive)
    meta=check['full_pole_and_snapshot']
    if 'scalar_div' in meta or meta['grid']!=grid:
        raise ArithmeticError('Unexpected refined scalar grid or mismatched band')
    G=[geometry(p) for p in pair]
    for mesh,(e,_) in enumerate(G,1):
        expected=topology(card,mesh,None,topology_only=True,mirror_azimuth=meta.get('mirror_azimuth',False))
        if not np.array_equal(e,expected['edges']) or any(expected[k]!=v for k,v in meta.items() if k!='mesh'):
            raise ArithmeticError('Native support or pole metadata fails exact reconstruction')
    return G,sha,check


def run(fine,third,states,archive):
    x,y=[audit(p,s,archive,g) for p,s,g in zip((fine,third),states,('fine','third'),strict=True)]
    result=[]
    for ((e,A),(f,B)) in zip(x[0],y[0],strict=True):
        r={}
        for channel,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
            r[channel]=compare(e,f,A[:,i,j],B[:,i,j])
        result.append(r)
    a,b=result
    defect=max(float(np.max(abs(np.asarray(a[k]['cdf_enclosures_at_cuts'])-t['cdf_enclosures_at_cuts']))) for k,t in b.items())
    ratio_defect=max(abs(a[k]['upper_over_lower_band_response']/t['upper_over_lower_band_response']-1.) for k,t in b.items())
    sharp=max(float(np.max(abs(np.asarray(a[k]['sharp_compatible_TV']['minimum_TV_bounds'])-
                              t['sharp_compatible_TV']['minimum_TV_bounds']))) for k,t in b.items())
    if max(defect,ratio_defect)>6e-5 or sharp>2e-9:
        raise ArithmeticError('Two q meshes disagree on native CDF enclosures or response ratio')
    return dict(x=.05,charge_subcells=2,grids=['fine','third'],steps=2048,
                states_sha256_uncompressed=[x[1],y[1]],gates=[x[2],y[2]],
                exact_shared_q_edges=np.intersect1d(x[0][1][0],y[0][1][0]).tolist(),
                maximum_paired_cdf_defect=defect,maximum_relative_paired_response_ratio_defect=ratio_defect,
                maximum_absolute_paired_sharp_tv_defect=sharp,
                observables=b,scope='Non-nested held finite band profiles, no within-panel interpolation or physical convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('fine','third','states'):
        p.add_argument('--'+name,nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive),indent=2,allow_nan=False))
