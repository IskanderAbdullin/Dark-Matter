"""Symmetric network splitting on individually exact number-conserving Bose edges."""
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from finite_mass_qed_bose_cell_exchange_reference import exchange,rates


def entropy(x):
    return np.log1p(x)+x*np.log1p(1/x) if x>0 else 0.


def split(y,c,A,B,energy,t,m,compensate=True):
    y0=np.asarray(y,dtype=float)
    u=y0.copy()
    initial=np.dot(c,u)
    parts=[[] for _ in A]
    diss=[]; scalar_extent=[]
    for step in range(m):
        for j in list(range(len(A)))+list(range(len(A)-1,-1,-1)):
            n=j+1
            prev=u[[0,n]].copy()
            u[0],u[n],_=exchange(*prev,c[0],c[n],A[j],B[j],t/(2*m))
            extent=c[n]*(prev[1]-u[n])
            parts[j].append(extent)
            scalar_extent.append(extent)
            if compensate:
                # Reconstruct from the original high part plus the fsum
                # of raw events. Keep the signed event ledger separately.
                u[0]=fsum((c[0]*y0[0],fsum(scalar_extent)))/c[0]
            exact_entropy=(c[0]*(entropy(u[0])-entropy(prev[0]))+
                           c[n]*(entropy(u[n])-entropy(prev[1]))+
                           extent*np.log(A[j]/B[j]))
            diss.append(exact_entropy)
            if max(abs(np.dot(c,u)/initial-1.),float(-u.min()))>8e-14:
                raise ArithmeticError('Perturbed global count or negative population')
    e=[fsum(q) for q in parts]
    # Stored reaction extent, not the difference of two nearly equal scalar states.
    E=np.r_[fsum(e)*energy[0],-np.array(e)*energy[1:],
             fsum(d*(e0-energy[0]) for d,e0 in zip(e,energy[1:],strict=True))]
    return u,np.array(e),fsum(E),min(diss)


def dilute():
    y=np.array([.45,.9,.5,.3,1.2])
    c=np.array([.8,1.2,1.5,.3,2.])
    A=np.array([.23,.42,.31,.59])
    B=np.array([.47,.15,.81,.42])
    Q=np.array([.65,.71,.83,.91,1.33])
    T=1.2
    def rhs(_,z):
        rr=np.array([rates(z[0],z[j+1],A[j],B[j]) for j in range(len(A))])
        return np.r_[sum(rr)/c[0],-rr/c[1:],rr]
    ref=solve_ivp(rhs,[0,T],np.r_[y,np.zeros(len(A))],rtol=9e-14,atol=2e-15,method='DOP853').y[:,-1]
    errors=[]
    for m in (8,16,32,64):
        x,I,e,s=split(y,c,A,B,Q,T,m)
        errors.append(float(max(abs(np.r_[x,I]-ref))))
        if abs(e)>2e-16 or s< -1e-13:
            raise ArithmeticError(('Split heat or entropy',m,e,s))
    ratio=[a/b for a,b in zip(errors[:-1],errors[1:],strict=True)]
    if min(ratio)<3.7:raise ArithmeticError(('Lost 2nd order',ratio))
    Y=1/np.expm1(Q/.83)
    B0=A*np.exp((Q[0]-Q[1:])/.83)
    X,I,E,S=split(Y,c,A,B0,Q,.8,8)
    if max(abs(X/Y-1.))>2e-13 or max(abs(I))>2e-13:
        raise ArithmeticError('A multichannel shared thermal fixed point moved')
    return dict(meshes=[8,16,32,64],absolute_errors=errors,convergence_ratios=ratio,
                thermal_state_defect=max(abs(X/Y-1.)))


def archived():
    root=Path(__file__).parent/'results'
    p=json.loads((root/'finite-mass-qed-timelike-pole-affine-source-2026-11-28.json').read_text())['data']
    out=[]
    import gzip
    for g in ('fine','third'):
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as s:
            d=json.load(s)
        q=4.1e-6
        n=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],q)-1]
        E=np.hypot(d['checkpoints'][-1]['scalar_pole_over_T'],q)
        v=[r for r in p if r['grid']==g]
        A,B,_=np.array([r['raw_kernel'] for r in v]).T
        y=np.array([n,.08,.15,.4,.3])
        C=np.array([1.,.2,1.,1.2,3.])
        H=np.r_[E,[r['omega'] for r in v]]
        t=5/(n*max(A))
        x,d,e,s=split(y,C,A,B,H,t,24)
        v0=fsum(d)
        uncontrolled=split(y,C,A,B,H,t,24,False)[0]
        if abs(e)>2e-17 or s< -1e-11:raise ArithmeticError((g,e,s))
        out.append(dict(grid=g,total_scalar_extent=v0,native_compensated_rounding=x[0]-y[0],
                        naive_repeated_scalar_update=uncontrolled[0]-y[0],
                        photon_final=x[1:].tolist(),energy_roundoff=abs(e),
                        compensated_charge_energy=fsum(I*(w-E) for I,w in zip(d,H[1:],strict=True))))
    return out


if __name__=='__main__':
    print(json.dumps(dict(scope='Multiple held reservoir edges, not a self-consistent relaxation model.',
                          controls=dilute(),archived_shapes=archived()),indent=2,allow_nan=False))
