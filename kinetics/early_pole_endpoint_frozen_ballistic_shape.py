"""Isotropic constant-velocity characteristic of positive frozen scalar archives."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_frozen_stress_sweep import weights,nodal_frozen
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def loss(z):
    """Stable angular average of 1-cos(z*mu), mu uniform on [-1,1]."""
    z=np.asarray(z,dtype=float)
    q=z*z
    series=q*(1/6+q*(-1/120+q*(1/5040+q*(-1/362880))))
    return np.where(abs(z)<=.045,series,1-np.sinc(z/np.pi))


def measure(q,w,m):
    if m<0 or np.any(q<=0) or min(w)<0 or not np.isfinite(m) or not np.all(np.isfinite([q,w])):
        raise ArithmeticError('Need finite positive-velocity measured archive')
    v=q/np.hypot(q,m)
    n=fsum(w)
    if n<=0:
        raise ArithmeticError('Empty archive has no characteristic function')
    p=w/n
    a,b=(fsum(p*v**k) for k in (2,4))
    if not 0<a<=1+1e-15 or not a*a<=b+2e-15 or b>a+1e-15:
        raise ArithmeticError('Velocity variance and fourth moment fail their probability bounds')
    return v,p,a,b


def eval_meas(v,p,phase):
    g=loss(v*phase)
    ans=fsum(p*g)
    upper=phase*phase*fsum(p*v*v)/6
    if not 0<=ans<=min(upper+2e-15,2.):
        raise ArithmeticError('Angular probability loss outside bounds')
    return ans


def run(native,geometric,gauss):
    rows=[unpack(c,bool(i),.05) for i,c in enumerate((native,*geometric))]
    d0,z0=rows[0]
    M=d0.state(.05,z0)[0]
    charts=[]
    for j,(d,y) in enumerate(rows):
        tag=1 if j==0 else geometric[j-1]['warm_div']
        q=[weights(d,y,n) for n in (48,96,192)]
        charts.append((f'W{tag}',q))
    for state in gauss:
        rep,f,k=nodal_frozen(state,d0,M)
        q=[weights(d0,f,n) for n in (48,96,192)]
        if abs(fsum(q[-1][1])/rep[0]-1)>3e-14:
            raise ArithmeticError('Lost archive scalar count')
        charts.append((f'G{k}',q))
    labels=[n for n,_ in charts]
    if labels!=['W1','W4','G2','G3','G4']:
        raise ArithmeticError('Use independent W1/W4/G2/G3/G4 cases')
    output=[]
    for m in (0.,float(M),.01,1.,10.,100.):
        a={name:[measure(q,w,m) for q,w,_ in triplet] for name,triplet in charts}
        lag=max(abs(v[2]/w[2]-1) for p in a.values() for v,w in zip(p[:-1],p[1:],strict=True))
        if lag>2e-10:
            raise ArithmeticError('Unresolved second velocity moment')
        xp=[.001,.003,.01,.03,.1,.3,1.,2.]
        orders=[]
        reject=[]
        for u in xp:
            h={}
            for test in ('same_native_phase','matched_rms_phase'):
                arr=[]
                for k,row in enumerate(a.values()):
                    L=[u/np.sqrt((a['W1'][i] if test=='same_native_phase' else row[i])[2]) for i in range(len(row))]
                    A=[eval_meas(v,p,s) for (v,p,_,_),s in zip(row,L,strict=True)]
                    ds=abs(A[-2]-A[-1])
                    if m==0 and abs(A[-1]-float(loss(u)))>4e-15:
                        raise ArithmeticError('Massless angular law did not reproduce the exact answer')
                    V,P,*_=row[-1]
                    mask=charts[k][1][-1][2]
                    warm=eval_meas(V[mask],P[mask],L[-1])
                    other=eval_meas(V[~mask],P[~mask],L[-1])
                    if abs(warm+other-A[-1])>3e-15:
                        raise ArithmeticError('Warm/cold partition lost angular response')
                    arr.append((A[-1],ds,L[-1],warm,other))
                if max(x[1] for x in arr)>max(3e-13,3e-9*max(x[0] for x in arr)):
                    reject.append(dict(u=u,phase_convention=test,max_absolute_96_192_defect=max(x[1] for x in arr)))
                    continue
                value,paired,L,W,C=map(np.asarray,zip(*arr,strict=True))
                h[test]=dict(loss=value.tolist(),raw_phase=L.tolist(),loss_difference_from_W1=(value-value[0]).tolist(),
                             partial_loss_q_ge_0p01=W.tolist(),partial_loss_q_lt_0p01=C.tolist(),
                             G4_minus_G3=float(value[-1]-value[-2]),paired_96_192_defect=paired.tolist())
            if h:
                orders.append(dict(u=u,**h))
        series=[]
        for label,row in a.items():
            v,p,second,fourth=row[-1]
            scale=.025/max(v)
            x=eval_meas(v,p,scale)
            predicted=scale**2*second/6-scale**4*fourth/120
            lim=scale**6*fsum(p*v**6)/5040
            if not -3e-16<=x-predicted<=lim+3e-16:
                raise ArithmeticError('Ballistic moment expansion failed')
            series.append(dict(chart=label,second=second,fourth=fourth,
                               kurtosis_velocity_squared=fourth/second**2,analytic_small_phase_check_abs=abs(x-predicted)))
        output.append(dict(formal_mass=m,measured_velocity_moments=series,accepted=orders,
                           unresolved=reject,largest_rel_second_moment_mesh_defect=lag))
    return dict(x=.05,chart_order=labels,measurement_rows=output,
                scope='Frozen positive spectra and formal observation mass; collisionless constant velocities, no expansion, sources, transfer function or physical k mapping.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--w4',type=Path,required=True)
    p.add_argument('--gauss',nargs=3,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),[read(a.w4)],[read(f) for f in a.gauss]),indent=2,allow_nan=False))
