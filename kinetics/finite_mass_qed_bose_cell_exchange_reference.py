"""Exact conditional two-mode bosonic exchange with stable minority occupation."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def root(a,b,c,D):
    p=c-a-b
    return (p+D)/(2*c) if p>0 else 2*a/(D-p)


def one_plus_factor(a,b,c,D,x0):
    if c==0:
        return 1.
    p=c-a-b
    at_zero=(2*a*c)/(D*(D+p)) if p>0 else (D-p)/(2*D)
    return at_zero+c*x0/D


def exchange(ns,np0,vs,vp,A,B,t):
    """Hold charged bath coefficients; evolve the Bose pair for a single edge."""
    if not np.all(np.isfinite([ns,np0,t,vs,vp,A,B])) or min(ns,np0,t)<0 or min(vs,vp,A,B)<=0:
        raise ValueError('Unphysical population, capacity, timespan or direction')
    Ns,Np=vs*ns,vp*np0
    total=Ns+Np
    a,b=A/vp,B/vs
    c=(A-B)*total/(vs*vp)
    D=np.hypot(c+a-b,2*np.sqrt(a)*np.sqrt(b))
    if t==0 or total==0:
        return ns,np0,D
    x,y=Ns/total,Np/total
    rx,ry=root(a,b,c,D),root(b,a,-c,D)
    rt=rx+ry
    rx,ry=rx/rt,ry/rt
    H=one_plus_factor(a,b,c,D,x) if c>=0 else one_plus_factor(b,a,-c,D,y)
    e=np.exp(-D*t)
    r=-np.expm1(-D*t)
    den=e+r*H
    sx=(e*x+r*H*rx)/den
    px=(e*y+r*H*ry)/den
    s=(sx+px)
    return (sx/s)*(total/vs),(px/s)*(total/vp),D


def log_affinity(S,P,A,B):
    return np.log(A)-np.log(B)+np.log(P)-np.log1p(P)-np.log(S)+np.log1p(S)


def rates(S,P,A,B):
    return A*P*(1+S)-B*S*(1+P)


def fairly():
    errs=[];balances=[];entropy=[];split=[]
    for A,B,vs,vp in ((1.,1.,1.,1.),(.1,5.,2.,.25),(6.,.07,.25,3.),(.4,.7,.5,1.6)):
        for S,P in ((.1,1.4),(.3,.02),(2.,.3),(0.,.8),(.1,0.)):
            for t in (.00002,.17,1.4):
                def rhs(_,y):
                    e=rates(*y,A,B)
                    return [e/vs,-e/vp]
                ref=solve_ivp(rhs,[0,t],[S,P],method='DOP853',rtol=3e-13,atol=2e-15)
                s,p,D=exchange(S,P,vs,vp,A,B,t)
                errs.append(max(abs(np.array([s,p])-ref.y[:,-1])))
                balances.append(abs((vs*s+vp*p)/(vs*S+vp*P)-1.))
                if S and P:
                    l0,l1=log_affinity(S,P,A,B),log_affinity(s,p,A,B)
                    if l1*l0< -5e-13 or abs(l1)>abs(l0)+1e-13:
                        raise ArithmeticError('Exact isolated-edge flow crossed equilibrium')
                    entropy.append(rates(s,p,A,B)*l1)
                a,b,_=exchange(S,P,vs,vp,A,B,t*.4)
                z,y,_=exchange(a,b,vs,vp,A,B,t*.6)
                split.append(max(abs(z-s),abs(y-p)))
    if max(errs)>3e-12 or max(balances)>8e-15 or min(entropy)<-1e-15 or max(split)>8e-15:
        raise ArithmeticError((max(errs),max(balances),min(entropy),max(split)))
    return dict(cases=len(errs),direct_dop853_max_abs_defect=max(errs),pair_number_max_rel_defect=max(balances),
                semigroup_max_abs_defect=max(split),min_instant_entropy_production=min(entropy))


def archive():
    root=Path(__file__).parent/'results'
    p=json.loads((root/'finite-mass-qed-timelike-pole-affine-source-2026-11-28.json').read_text())['data']
    out=[]
    for z in p:
        g=z['grid']
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:d=json.load(f)
        S=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],4.1e-6)-1]
        A,B,_=z['raw_kernel']
        for ratio in (.01,1.,100.):
            Pg=.3
            _,_,D=exchange(S,Pg,1.,ratio,A,B,1e-30)
            s,p,D=exchange(S,Pg,1.,ratio,A,B,60/D)
            ell=log_affinity(s,p,A,B)
            if abs(ell)>2e-11 or min(p,s)<=0. or abs((s+ratio*p)/(S+ratio*Pg)-1)>4e-15:
                raise ArithmeticError(('Large occupation minority lost',g,z['polarization'],z['k'],ratio,ell,p))
            es=np.hypot(d['checkpoints'][-1]['scalar_pole_over_T'],4.1e-6)
            w=z['omega']
            # Once the raw scalar occupation is near 1e16, subtracting two
            # 64-bit scalar populations need not return the reaction extent.
            extent=ratio*(Pg-p)
            from math import fsum
            energy=[es*extent,-w*extent,(w-es)*extent]
            residual=abs(fsum(energy))
            if residual>3e-16*max(abs(np.asarray(energy))):
                raise ArithmeticError('Reaction energy bookkeeping did not close')
            out.append(dict(grid=g,k=z['k'],polarization=z['polarization'],capacity_ratio=ratio,
                            minority_photon_occupation=p,late_log_affinity_defect=abs(ell),
                            scaled_deadline=60.,same_fixed_scalar_mode=4.1e-6,
                            extent_from_photon=extent,
                            direct_scalar_update_error=(s-S)-extent,energy_ledger=energy,energy_defect=residual))
    return out


def common_heat():
    er=[]
    for T in (.3,.9,1.5):
        for es,w in ((.00087,.0023),(.001,.0018),(.04,.15)):
            nS,nP=(1/np.expm1(E/T) for E in (es,w))
            B=np.exp((es-w)/T)
            s,p,_=exchange(nS,nP,.8,1.7,1.,B,90.)
            er.extend((abs(s/nS-1),abs(p/nP-1)))
    if max(er)>4e-12:
        raise ArithmeticError('Common Bose fixed point in the held model')
    return max(er)


if __name__=='__main__':
    print(json.dumps(dict(scope='Single reversible edge with held charged reservoirs, not a population model.',
                         independent_tests=fairly(),thermal_fixed_point_max_rel=common_heat(),archived_shapes=archive()),indent=2,allow_nan=False))
