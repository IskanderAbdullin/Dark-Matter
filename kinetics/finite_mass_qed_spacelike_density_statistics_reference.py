"""Free charge-density statistics with a positive recoil integral and a Ward gate."""

import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit, roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import norm_spinors, slash

X, A = roots_legendre(3)
GT = slash([1., 0., 0., 0.])
GZ = slash([0., 0., 0., -1.])


def aperture(m,k,w):
    a = (k-w)*(k+w)/(k*k)
    s = a*k*k
    h = 2*m/np.sqrt(a)
    b = np.hypot(k,h)
    return a,s,b,(s+h*h)/(2*(b+w))


def weight(E,m,p):
    a,_,b,e0 = p
    return 8*(m*m/a+(E-e0)*(E-e0+b))


def exact(m,e,f,k,w):
    p = aperture(m,k,w)
    a,_,b,e0 = p
    if e[-1] <= e0:
        return np.zeros(4)
    g = np.unique(np.clip(np.r_[e0,e[-1],e,e-w],e0,e[-1]))
    q,mid = np.diff(g)/2,(g[:-1]+g[1:])/2
    weights = q*(weight(mid[:,None]+q[:,None]*X,m,p) @ A)
    ext = np.r_[0.,f,0.]
    l,r = (ext[np.searchsorted(e,t)] for t in (mid,mid+w))
    G,L,D = weights @ (l*(1-r)), weights @ (r*(1-l)), weights @ (l-r)
    x,y = (np.clip((t-e0)/w,0,1) for t in (e[:-1],e[1:]))
    clip = np.minimum(np.maximum(e[:-1],e0+w),e[1:])
    low = 8*(y-x)*(m*m/a+b*w*(x+y)/2+w*w*(x*x+x*y+y*y)/3)
    stable = f @ (low+8*(e[1:]-clip)*(e[1:]+clip))
    return np.array([G,L,D,stable])


def soft(m,e,f,u):
    p = m/np.sqrt((1-u)*(1+u))
    H = e[1:]
    L = np.minimum(np.maximum(e[:-1],p),H)
    N = (f*(1-f)) @ ((H-L)*(H*H+H*L+L*L))/3
    ext = np.r_[0.,f,0.]
    left,right = (ext[np.searchsorted(e,p,s)] for s in ('left','right'))
    D = (left+right)*p*p/2 + f @ ((H-L)*(H+L))
    jumps = (ext[:-1]-ext[1:])*np.where(e>p,1.,np.where(e==p,.5,0.))
    return np.array([N,D,jumps @ (e*e)])


def spin_gate(m,k,w):
    p=aperture(m,k,w)
    _,s,_,e0=p
    u=w/k
    out=0.
    for E in (e0+.03,e0+.1,e0+1.):
        pz=(w*E-s/2)/k
        px=np.sqrt((E-m)*(E+m)-pz*pz)
        P=np.array([E,px,0.,pz])
        B=P+np.array([w,0.,0.,k])
        tt,zz,plus,minus,ward=(norm_spinors(P,B,m,m,g) for g in (GT,GZ,GT+GZ,GT-GZ,u*GT-GZ))
        out=max(out,abs(tt/weight(E,m,p)-1),abs(zz/tt-u*u),abs((plus-minus)/(4*tt)-u),ward/tt)
    return out


def run():
    data=[]; ward=[]; thermal=[]
    for grid in ('fine','third'):
        root=Path(__file__).parent/'results'
        with gzip.open(root/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as stream:
            d=json.load(stream)
        m=d['x_end']
        c=[(q,np.hypot(q,d[f'endpoint_{n}_momentum_edges']),np.asarray(d[f'endpoint_{n}_occupations']))
           for q,n in ((m,'parent'),(m*d['mass_ratio'],'fermion'))]
        rows=[]; check=[]
        for u in (.0047,.25,.5,.9,.99,.9996):
            S,D,J=np.sum([soft(*q,u) for q in c],axis=0)
            k=1e-12
            P=sum(exact(*q,k,k*u)[3] for q in c)
            check.append(max(abs(P/(8*D)-1.),abs(J/D-1.)))
            for k in (.0002,.0008,.00256,.01,.1):
                w=k*u
                G,L,R,P=np.sum([exact(*q,k,w) for q in c],axis=0)
                if min(G,L,R,P)<=0 or abs(R/(w*P)-1)>2e-8:
                    raise ArithmeticError((grid,k,u,G,L,R,P))
                ratio=(G+L)/(2*P)
                rows.append(dict(u=u,k=k,soft=S/D,exact=ratio,relative_recoil=ratio/(S/D)-1,
                                 stable_vs_direct=abs(R/(w*P)-1)))
                ward.extend(spin_gate(q[0],k,w) for q in c)
        data.append(dict(grid=grid,soft_limit_and_jump_max_defect=max(check),rows=rows))
    for m in (.003,.05):
        for k,u in ((.0002,.25),(.00256,.9),(.1,.99)):
            w=k*u; p=aperture(m,k,w); e0=p[-1]
            for T in (.4,1.):
                I=lambda h:quad(h,e0,e0+60*T)[0]
                G=I(lambda E:weight(E,m,p)*expit(-E/T)*expit((E+w)/T))
                L=I(lambda E:weight(E,m,p)*expit(-(E+w)/T)*expit(E/T))
                D=quad(lambda x:weight(e0+w*x,m,p)*expit(-(e0+w*x)/T),0,1)[0]
                D+=16*quad(lambda E:E*expit(-E/T),e0+w,e0+60*T)[0]
                thermal.append(max(abs(L/G/np.exp(-w/T)-1.),abs((G+L)/D*np.tanh(w/(2*T))/w-1.)))
    # Put a genuine nonmonotone jump exactly at the soft angular threshold.
    m,u,k=.25,.75,1e-7
    e0=m/np.sqrt((1-u)*(1+u))
    e,f=np.array([m,e0,2*e0,4*e0]),np.array([.2,.85,.1])
    _,D,J=soft(m,e,f,u)
    edge=max(abs(J/D-1.),abs(exact(m,e,f,k,k*u)[3]/(8*D)-1.))
    if edge>2e-8:
        raise ArithmeticError(('Threshold jump',edge))
    if max(ward)>2e-12 or max(thermal)>2e-12 or max(r['soft_limit_and_jump_max_defect'] for r in data)>3e-10:
        raise ArithmeticError((max(ward),max(thermal),[r['soft_limit_and_jump_max_defect'] for r in data]))
    return dict(scope='Free spacelike charge-density cut; not a longitudinal collective population.',
                density_spinor_ward_max_rel=max(ward),equilibrium_max_rel=max(thermal),
                planted_threshold_jump_max_rel=edge,data=data)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
