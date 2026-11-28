"""Compare exact G8/G9/mixed300 mature endpoints under prescribed free streaming."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_frozen_ballistic_shape import eval_meas
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_radiation_clock_streaming_probe import average_speed,clock_check
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import fit,scaled_energy,shell_rule
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def populations(card,base,m,control):
    tag,d,k=chart_for(card,budget_control=control)
    held=gauss(card,base,m,.05,expected_steps=2048,budget_control=control)[0]
    endpoints=d.ds.ed
    y=np.asarray(card['nodal_occupation'])
    N=y[:d.ns]*d.ds.w
    e=np.r_[0,np.cumsum(k)]
    prepared=[]
    for j,(a,b) in enumerate(zip(endpoints[:-1],endpoints[1:],strict=True)):
        q=d.ds.p[e[j]:e[j+1]]
        w=N[e[j]:e[j+1]]
        mass=fsum(w)
        if not mass:
            continue
        if k[j]==1:
            prepared.append((a,b,q,w,mass,None))
        else:
            u=fsum(w*scaled_energy(q,a,b,m))/mass
            prepared.append((a,b,q,w,mass,fit(a,b,m,mass,u)))
    reps=[]
    for n in (64,128):
        modes=[([],[]),([],[])]
        for a,b,q,w,z,F in prepared:
            Q,v,L=shell_rule(a,b,m,n)
            if F is None:
                out=np.exp(L)
                out=out*(z/fsum(out))
                branch=((Q,out),(Q,out))
            else:
                out=np.exp(L+F['log_occupation_at_lower_edge']+F['tilt']*v)
                if abs(fsum(out)/z-1.)>5e-12 or abs(fsum(out*v)/z-F['scaled_mean_energy'])>5e-12:
                    raise ArithmeticError('Retained warm count and energy did not replay')
                branch=((q,w),(Q,out))
            for (X,W),(x,v) in zip(modes,branch,strict=True):
                X.extend(x); W.extend(v)
        pack=[]
        for q,w in modes:
            X,W=np.array(q),np.array(w)
            total=fsum(W)
            if abs(total/held[0]-1.)>6e-14:
                raise ArithmeticError('Frozen streaming representation lost accepted scalar count')
            pack.append((X,W/total))
        reps.append(pack)
    return tag,held,reps


def triplet(pops,m,A,u,level,mode):
    data=[]
    for _,_,p in pops:
        q,w=p[level][mode]
        data.append((average_speed(q,m,A),w))
    s=np.array([fsum(v*v*w) for v,w in data])
    b=u/np.sqrt(s)
    R=np.array([[eval_meas(v,w,t) for t in b] for v,w in data])
    out=[]
    for j,k in ((0,1),(0,2),(1,2)):
        d=R[k,k]-R[j,j]
        fixed=(R[k,j]-R[j,j]+R[k,k]-R[j,k])/2
        phase=(R[k,k]-R[k,j]+R[j,k]-R[j,j])/2
        if abs(d-fixed-phase)>3e-15:
            raise ArithmeticError('Cannot decompose equal-rms contrast')
        out.append([R[k,0]-R[j,0],d,fixed,phase])
    return s,b,np.array(out),R


def run(native,accepted,control,state):
    (a,ha)=read(native)
    raw=accepted.read_bytes()
    lock=json.loads(raw)
    card=[read(f) for f in state]
    r=lock['rows'][-1]
    if r['x']!=.05 or lock['charts']!=['G8','G9','mixed300'] or r['states_sha256']!=[h for _,h in card]:
        raise ArithmeticError('Three exact stored x=.05 controls required')
    d,y=unpack(a,False,.05)
    m=d.state(.05,y)[0]
    if m!=r['common_native_mass']:
        raise ArithmeticError('No new mass fit is permitted')
    P=[populations(c,d,m,f) for (c,_),f in zip(card,(None,None,control),strict=True)]
    if [p[0] for p in P]!=[8,9,'mixed300']:
        raise ArithmeticError('Wrong chart ordering')
    v=np.array([p[1] for p in P])
    if np.max(abs(v[1:]/v[0]-1.-r['ratios_to_G8']))>2e-15:
        raise ArithmeticError('Existing held pressure comparison did not reproduce')
    rows=[]
    for A in (1.,2.,10.,100.,1000.):
        checksum=max([clock_check(p[2][-1][i][0][::83],m,A) for p in P for i in (0,1)] if A>1. else [0.])
        if checksum>5e-12:
            raise ArithmeticError('Analytic and independent integrated distance differ')
        out=[]
        for u in (.003,.03,.3,1.,2.):
            comparisons=[]
            for k in (0,1):
                lo,hi=[triplet(P,m,A,u,n,k) for n in (0,1)]
                s,b,T,R=hi
                corr=max(np.max(abs(lo[i]-hi[i])) for i in (2,3))
                if corr>max(3e-13,3e-9*np.max(R)):
                    raise ArithmeticError(('Independent frozen quadrature did not resolve',A,u,k,corr))
                comparisons.append(dict(mean_square_speed=s.tolist(),phase=b.tolist(),matrix=R.tolist(),
                                        differences=T.tolist(),maximum_paired_rule_difference=corr))
            out.append(dict(u=u,representations=comparisons,
                            second_minus_first_differences=(np.array(comparisons[1]['differences'])-comparisons[0]['differences']).tolist()))
        rows.append(dict(A=A,maximum_independent_clock_defect=checksum,rows=out))
    return dict(x=.05,formal_initial_mass=m,charts=lock['charts'],pair_order=[['G9','G8'],['mixed300','G8'],['mixed300','G9']],
                representations=['warm nodal atoms','warm positive count-and-energy fit'],
                difference_fields=['common_G8_raw_phase','equal_rms','symmetric_fixed_phase','symmetric_phase_rescaling'],
                rows=rows,input_sha256=[ha,hashlib.sha256(raw).hexdigest(),hashlib.sha256(control.read_bytes()).hexdigest(),*[h for _,h in card]],
                scope='Frozen characteristics in an imposed clock with no metric forcing or kinetic evolution. Two admissible warm reconstructions only, not a bound, spatial convergence, or a physical transfer function.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','accepted','control'):
        p.add_argument(k,type=Path)
    p.add_argument('state',type=Path,nargs=3)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
