"""Keep shell number and energy while smoothing the bare virtual overlap."""
import gzip,json
from functools import lru_cache
from math import fsum
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_closed_radial_pv_reference import angular_kernel, full_radial
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes

D=np.longdouble
X=Polynomial([0.,1.])
bump=X**2*(1-X)**2*(2*X-1)
GL=lru_cache(None)(roots_legendre)


class Smoother:
    def __init__(self,edges,occ,mass,width):
        self.m=D(mass)
        self.u=np.asarray(edges,dtype=D)**3/3
        self.mu=np.asarray(occ,dtype=D)
        x,w=GL(40)
        z,ww=(x+1)/2,w/2
        def energy(u):
            return np.hypot(np.cbrt(3*u),self.m)
        self.k=[]; self.edges=[]; self.used=[]; self.dual_defect=[]; self.extrema=[]
        edge=[occ[0]]
        for a,b in zip(occ[:-1],occ[1:],strict=True):
            if a in (0,1) and b in (0,1) and a!=b:
                raise ArithmeticError('Incompatible touching saturated cells.')
            edge.append(a if a in (0,1) else b if b in (0,1) else (a+b)/2)
        edge.append(occ[-1])
        for j, F in enumerate(self.mu):
            l,r=self.u[j:j+2]
            v,t=map(D,edge[j:j+2])
            cap=D(width)
            if F in (0,1):
                self.extrema.append((float(F),float(F)))
                fmid=F;c=0.
                cap=D(0)
            else:
                for _ in range(19):
                    fmid=(F-cap*(v+t)/2)/(1-cap)
                    def base(q):
                        n=np.clip(q/cap,0.,1.)
                        p=np.clip((1-q)/cap,0.,1.)
                        A=1-n*n*(3-2*n)
                        B=1-p*p*(3-2*p)
                        return fmid+(v-fmid)*A+(t-fmid)*B
                    d=r-l
                    temp=energy(l+d/2)
                    def piece(n):
                        x,y=GL(n)
                        grid=[0.,float(cap),float(1-cap),1.]
                        ab=[]
                        for A,B in zip(grid[:-1],grid[1:],strict=True):
                            q=D(A)+D(B-A)*(np.asarray(x,dtype=D)+1)/2
                            ab.append((q,D(B-A)*np.asarray(y,dtype=D)/2))
                        return ab
                    def weights(n,c):
                        n0=n1=D(0)
                        for q,y in piece(n):
                            val=base(q)+c*bump(q)-F
                            n0+=np.dot(y,val)
                            n1+=np.dot(y,val*(energy(l+d*q)-temp))
                        return n0,n1
                    q=np.array(z,dtype=D)
                    den=np.dot(ww,bump(q)*(energy(l+d*q)-temp))
                    c=-weights(40,0.)[1]/den
                    # The three pieces are degree at most five. Check the full
                    # capacity by endpoints plus all interior derivative roots.
                    cubic=1-X*X*(3-2*X)
                    CL=float(cap)
                    ranges=[(0.,CL,float(v-fmid)),(CL,1-2*CL,0.),(1.,-CL,float(t-fmid))]
                    feasible=True
                    for offset,extent,transition in ranges:
                        site=offset+extent*X
                        curve=Polynomial([float(fmid)])+transition*cubic+float(c)*bump(site)
                        q=[0.,1.,*[rr.real for rr in curve.deriv().roots()
                                    if abs(rr.imag)<1e-10 and 0<rr.real<1]]
                        q=D(offset)+D(extent)*np.asarray(q,dtype=D)
                        vals=base(q)+c*bump(q)
                        if min(vals)<-2e-16 or max(vals)>1+2e-16:
                            feasible=False
                            break
                    if feasible:
                        e0,e1=weights(68,c)
                        e1=e1/max(energy(r)-energy(l),D('1e-30'))
                        if max(abs(e0),abs(e1))>4e-12:
                            raise ArithmeticError('Two-moment conservative reconstruction did not match independent integration.')
                        self.dual_defect.append([float(abs(e0)),float(abs(e1))])
                        tape=np.linspace(0,1,129,dtype=D)
                        q=np.r_[cap*tape,cap+(1-2*cap)*tape,1-cap*tape]
                        sampled=base(q)+c*bump(q)
                        if np.min(sampled)<-2e-16 or np.max(sampled)>1+2e-16:
                            raise ArithmeticError('A separate sampled capacity bound failed.')
                        self.extrema.append((float(min(sampled)),float(max(sampled))))
                        break
                    cap/=2
                else:
                    raise ArithmeticError('No admissible continuous conservative cell at this requested width.')
            self.k.append((cap,fmid,v,t,c))
            self.used.append(float(cap))
            self.edges.extend([l,l+(r-l)*cap,r-(r-l)*cap,r])
        self.boundaries=np.hypot(np.cbrt(3*np.asarray(sorted(set(self.edges)),dtype=D)),self.m)

    def value(self,E):
        if E<self.m: return D(0)
        p2=(E-self.m)*(E+self.m)
        u=p2*np.sqrt(p2)/3
        j=np.searchsorted(self.u,u,side='right')-1
        if j<0 or j>=len(self.mu): return D(0)
        cap,mid,l,r,C=self.k[j]
        if cap==0:return mid
        q=(u-self.u[j])/(self.u[j+1]-self.u[j])
        h=np.clip(q/cap,0.,1.);g=np.clip((1-q)/cap,0.,1.)
        # Direct bounded evaluation. No clipping of the reconstructed result.
        return mid+(l-mid)*(1-h*h*(3-2*h))+(r-mid)*(1-g*g*(3-2*g))+C*bump(q)


def calc(d,E,cut,lo,hi,width,step,prepared=None,orders=(37,55)):
    M=d['x_end'];m=d['mass_ratio']*M
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(D(s),D(Q))
    name,ma,num=('parent',M,D(E)) if cut=='parent' else ('fermion',m,D(E)-es)
    rec=prepared if prepared is not None else Smoother(d[f'endpoint_{name}_momentum_edges'],d[f'endpoint_{name}_occupations'],ma,width)
    fn=angular_kernel(M,m,s,E,cut)
    _,C,_=solve_axes(M,m,s,E,cut)
    coeff=D(C)*np.sqrt((num-ma)*(num+ma))*rec.value(num)/2
    coords=[float(q-num) for q in rec.boundaries]
    start,end=coords[0],coords[-1]
    extras=[v for v in (lo-E,hi-E,0.) if start<v<end]
    tests=sorted({*coords,*extras,*[v for dd in np.geomspace(1e-10,max(abs(start),end),45)
                                   for v in (dd,-dd) if start<v<end],
                      *[b+e for b in (lo-E,hi-E) for e in (-1e-6,-1e-7,1e-7,1e-6)
                        if start<b+e<end]})
    def integral(n):
        x,w=GL(n)
        p=(x+1)*np.pi/2
        A=np.sin(p/2)**2
        W=np.pi/4*w*np.sin(p)
        vals=[]
        for a,b in zip(tests[:-1],tests[1:],strict=True):
            if b<=a:continue
            pos=D(a)+(b-a)*A.astype(D)
            f=[]
            for dr in pos:
                e=num+dr
                occ=rec.value(e)
                r=np.sqrt(max(0.,(e-ma)*(e+ma)))
                f.append(float(fn(dr)*r*occ/2-coeff/dr))
            vals.append((b-a)*(W@f))
        return fsum(vals)+float(coeff*np.log(abs(end/start)))
    a,b=(integral(n) for n in orders)
    if abs(a-b)>2e-5:
        raise ArithmeticError(f'Doubly conservative virtual quadrature did not resolve at {width}: {a}, {b}.')
    return dict(requested_width=width,used_nonzero_min=min((x for x in rec.used if x>0),default=0.),
                reconstruction_shell_defects=np.asarray(rec.dual_defect).max(axis=0).tolist(),
                sampled_pauli_min_and_max=[min(x[0] for x in rec.extrema),max(x[1] for x in rec.extrema)],
                internal_continuous_occupation_at_coincidence=float(rec.value(num)),
                bare_virtual=b,radial_mesh_defect=abs(a-b))


def run():
    out=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=d['mass_ratio']*M;s=d['checkpoints'][-1]['scalar_pole_over_T']
        es,lo,hi,*_=intervals(M,s,m,Q)
        cuts=('parent','daughter')
        obs=[]
        for cut in cuts:
            E=(lo+hi)/2
            obs.append(dict(cut=cut,E=E,tuned_step=False,
                            bare_steps_and_linear=full_radial(d,E,cut,lo,hi,False),
                            conservative=[calc(d,E,cut,lo,hi,w,False) for w in (.015,.05,.12)]))
        cut='daughter' if grid=='fine' else 'parent'
        name,ma,shift=('fermion',m,es) if cut=='daughter' else ('parent',M,0.)
        edges=np.hypot(d[f'endpoint_{name}_momentum_edges'][1:-1],ma)+shift
        ok=np.flatnonzero((edges>lo+4e-4)&(edges<hi-4e-4))
        ix=ok[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[ok]))]
        E=float(edges[ix])
        obs.append(dict(cut=cut,E=E,tuned_step=True,
                        bare_steps_and_linear=full_radial(d,E,cut,lo,hi,True),
                        conservative=[calc(d,E,cut,lo,hi,w,True) for w in (.015,.05,.12)]))
        out.append(dict(grid=grid,cases=obs))
    return dict(scope='Continuous bounded occupations preserving the number and energy in each momentum cell; isolated bare positive virtual projection only, no physical rate.',grids=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
