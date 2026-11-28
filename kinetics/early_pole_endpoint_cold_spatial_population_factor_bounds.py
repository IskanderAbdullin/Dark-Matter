"""Bound explicit population substitution across misaligned cold scalar grids."""
import argparse
import json
from functools import reduce
from math import fsum,log
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry
from early_pole_endpoint_cold_shared_reference_scan import fit_refined_left


FIELDS=(('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0))


def crossed(q,bins,raw,source_edges,source_f,target_edges,target_f,inverse):
    """Positive-kernel panel bounds; unknown subdivisions are never assigned a rate."""
    lo,hi=q[:-1],q[1:]
    home=np.searchsorted(source_edges,lo,side='right')-1
    dest=np.searchsorted(bins,lo,side='right')-1
    begin=np.searchsorted(target_edges,lo,side='right')-1
    end=np.searchsorted(target_edges,hi,side='left')
    if (len(raw)!=len(lo) or np.any(raw<0) or not np.all(np.isfinite(raw)) or
            np.any(np.diff(q)<=0) or np.any(home<0) or np.any(home>=len(source_f))
            or np.any(dest<0) or np.any(dest>=len(bins)-1)
            or np.any(hi>source_edges[home+1]) or np.any(hi>bins[dest+1])
            or np.any(begin<0) or np.any(end>len(target_f)) or np.any(begin>=end)):
        raise ArithmeticError('Intervals do not bind to native panels and display bins')
    ft,fs=np.asarray(target_f),np.asarray(source_f)
    if not np.all(np.isfinite(ft)) or not np.all(np.isfinite(fs)) or min(ft)<0. or min(fs)<0.:
        raise ArithmeticError('Need nonnegative finite scalar populations')
    offset=0. if inverse else 1.
    factors=ft+offset
    if min(fs[home]+offset)<=0.:
        raise ArithmeticError('Cannot divide out a zero native explicit factor')
    K=raw/(fs[home]+offset)
    lim=np.array([[min(factors[a:b]),max(factors[a:b])] for a,b in zip(begin,end,strict=True)])
    vals=K[:,None]*lim
    out=np.stack([np.bincount(dest,weights=vals[:,i],minlength=len(bins)-1) for i in (0,1)],axis=1)
    return out,int(np.count_nonzero(lim[:,0]!=lim[:,1]))


def normalized_bounds(c):
    l,h=c.T
    if np.any(l<0) or np.any(h<l) or fsum(l)<=0.:
        raise ArithmeticError('Need positive independent bin bounds')
    return np.stack([l/(fsum(h)-h+l),h/(fsum(l)-l+h)],axis=1)


def tv_with_interval_profile(native,bound):
    hit=np.maximum(0.,native-bound[:,1])
    miss=np.maximum(0.,bound[:,0]-native)
    upper=fsum(np.maximum(abs(native-bound[:,0]),abs(native-bound[:,1])))/2
    guard=64*np.finfo(float).eps
    return [max(0.,fsum(hit)-guard,fsum(miss)-guard),min(1.,upper+guard)]


def field(books,cards,bins,i,j):
    B=np.zeros((2,2,len(bins)-1,2))
    counts=np.zeros((2,2),dtype=int)
    for s,(q,C) in enumerate(books):
        for t,(other,f,_,_) in enumerate(cards):
            origin,g,_,_=cards[s]
            B[s,t],counts[s,t]=crossed(q,bins,C[:,i,j],origin.ds.ed,g,other.ds.ed,f,j>0)
        ref=fsum(C[:,i,j])
        if counts[s,s] or ref<=0. or abs(fsum(B[s,s,:,0])/ref-1.)>4e-14:
            raise ArithmeticError('Native diagonal does not reconstruct exactly')
    sums=np.array([[[fsum(B[s,t,:,k]) for k in (0,1)] for t in (0,1)] for s in (0,1)])
    if np.min(sums)<=0.:
        raise ArithmeticError('Cannot bound log substitution with nonpositive gross responses')
    d=np.log(sums)
    total=d[1,1,0]-d[0,0,0]
    pop=np.array([(d[0,1,a]+d[1,1,0]-d[0,0,0]-d[1,0,1-a])/2. for a in (0,1)])
    p0,p1=[B[k,k,:,0]/sums[k,k,0] for k in (0,1)]
    z01,z10=normalized_bounds(B[0,1]),normalized_bounds(B[1,0])
    pop_change=[tv_with_interval_profile(p0,z01),tv_with_interval_profile(p1,z10)]
    kernel_change=[tv_with_interval_profile(p0,z10),tv_with_interval_profile(p1,z01)]
    exact={}
    if counts[1,0]==0:
        anchor=log(sums[1,0,0]/sums[0,0,0])
        update=log(sums[1,1,0]/sums[1,0,0])
        if abs(total-anchor-update)>2e-14:
            raise ArithmeticError('Population-anchored path does not reconstruct its endpoint contrast')
        exact=dict(exact_kernel_first_log_contrasts={'kernel_at_coarse_population':anchor,
                    'population_at_fine_kernel':update},
                   exact_kernel_first_population_share_of_native_log=(update/total if total else None))
    return dict(**exact,native_log_response_ratio=float(total),explicit_population_log_bounds=pop.tolist(),
                explicit_population_share_of_native_log_bounds=(sorted(pop/total) if total else None),
                combined_bare_kernel_log_bounds=(total-pop[::-1]).tolist(),
                native_profile_TV=fsum(abs(p1-p0))/2,
                same_kernel_population_profile_TV_bounds=pop_change,
                same_population_kernel_profile_TV_bounds=kernel_change,
                crossed_panels_with_nonconstant_factor=counts.tolist(),
                crossed_band_response_bounds=sums.tolist(),
                native_and_crossed_bin_response_bounds=B.tolist())


def run(pairs,states,archive,*,refined_s8=None,splits=None):
    levels=(2,4,8,16)
    cards=[fit_state(s,archive,n,r) for n,r,s in zip(levels,pairs,states,strict=True)]
    if (refined_s8 is None)!=(splits is None):
        raise ArithmeticError('Supply both refined S8 meshes and the exact plan')
    if refined_s8 is not None:
        check=fit_refined_left(refined_s8,pairs[2],states[2],archive,splits,level=8)
        if check[-1]!=cards[2][-1]:
            raise ArithmeticError('Refined panel S8 state changed')
    mesh=[]
    for t in (0,1):
        books=[geometry(p[t]) for p in pairs]
        bins=reduce(np.intersect1d,(p[0] for p in books))
        if refined_s8 is not None:
            books[2]=geometry(refined_s8[t])
        if [bins[0],bins[-1]]!=[1e-6,1.5e-5] or len(bins)<6:
            raise ArithmeticError('Incomplete exact common cold display')
        mesh.append((bins.tolist(),[dict(levels=list(levels[i:i+2]),observables={
            name:field(books[i:i+2],cards[i:i+2],bins,r,c) for name,r,c in FIELDS})
            for i in range(len(levels)-1)]))
    (a,x),(b,y)=mesh
    if a!=b:
        raise ArithmeticError('Exact common bins differed across q meshes')
    ratio_error=0.
    tv_error=0.
    for u,v in zip(x,y,strict=True):
        for key,L in u['observables'].items():
            R=v['observables'][key]
            A,B=[np.array(z['crossed_band_response_bounds']) for z in (L,R)]
            ratio_error=max(ratio_error,float(np.max(abs(A/B-1.))))
            for k in ('native_profile_TV','same_kernel_population_profile_TV_bounds',
                      'same_population_kernel_profile_TV_bounds','explicit_population_log_bounds'):
                tv_error=max(tv_error,float(np.max(abs(np.asarray(L[k])-R[k]))))
            if L['crossed_panels_with_nonconstant_factor']!=R['crossed_panels_with_nonconstant_factor']:
                raise ArithmeticError('Two meshes show different positive-factor geometry')
    if ratio_error>6e-5 or tv_error>2e-9:
        raise ArithmeticError('Paired substitution bounds did not resolve')
    mode={} if refined_s8 is None else dict(S8_target_population_edge_reintegration=True)
    return dict(**mode,x=.05,levels=list(levels),common_q_edges=b,adjacent=y,
                states_sha256_uncompressed=[c[-1] for c in cards],
                maximum_paired_fractional_response_bound_defect=ratio_error,
                maximum_paired_log_or_tv_bound_defect=tv_error,
                scope='Positive-panel outer bounds; no within-panel shape inferred, and bare kernel is not further decomposed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in (2,4,8,16):
        p.add_argument(f'--s{n}',type=Path,nargs=2,required=True)
    p.add_argument('--states',type=Path,nargs=4,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--refined-s8',type=Path,nargs=2)
    p.add_argument('--s8-split-plan',type=Path)
    arg=p.parse_args()
    splits=json.loads(arg.s8_split_plan.read_text())['targeted_parent_panel_splits'] if arg.s8_split_plan else None
    print(json.dumps(run((arg.s2,arg.s4,arg.s8,arg.s16),arg.states,arg.archive,refined_s8=arg.refined_s8,splits=splits),indent=2,allow_nan=False))
