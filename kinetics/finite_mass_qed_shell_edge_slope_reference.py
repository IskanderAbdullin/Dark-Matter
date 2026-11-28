"""The hard log pole derivative at a finite-volume occupation discontinuity."""
import functools
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from finite_mass_qed_archived_aperture_reference import Frozen, Q, solve_window


def build(edges,f,fraction,cut):
    """Return piecewise C1 smoothstep transitions in u=p³/3, or exact cell constants."""
    u=edges**3/3
    f=np.asarray(f)
    lo=cut**3/3
    seg=[]
    if fraction==0:
        for l,h,F in zip(u[:-1],u[1:],f,strict=True):
            if h>lo:
                seg.append((max(l,lo),h,l,h,F,F))
        return seg
    junction=np.r_[f[0],(f[:-1]+f[1:])/2,f[-1]]
    for j,(l,h,F) in enumerate(zip(u[:-1],u[1:],f,strict=True)):
        left,right=junction[j:j+2]
        if h<=lo:
            continue
        halfmean=(left+right)/2
        b=min(fraction,
              .5*F/halfmean if halfmean>0 else fraction,
              .5*(1-F)/(1-halfmean) if halfmean<1 else fraction)
        if b<=0:
            if left!=F or right!=F:
                raise ValueError('Saturated cell rejects these edge targets.')
            seg.append((max(l,lo),h,l,h,F,F))
            continue
        mid=(F-b*halfmean)/(1-b)
        w=(h-l)*b
        if not (0<=mid<=1):
            raise ValueError('Smooth conservative plateau violates capacity.')
        if abs((mid*(1-b)+b*halfmean)-F)>3e-16:
            raise ValueError('C1 reconstruction changed shell number.')
        for A,B,x,y in ((l,l+w,left,mid),(l+w,h-w,mid,mid),(h-w,h,mid,right)):
            if B>lo and B>A:
                seg.append((max(A,lo),B,A,B,x,y))
    return seg


class Pole:
    def __init__(self,edges,f,m,cut,fraction):
        self.m=m
        obj=np.asarray(build(edges,f,fraction,cut))
        self.l,self.r=np.cbrt(3*obj[:,:2]).T
        self.U,self.V,self.y,self.z=obj[:,2:].T
        self.yjump,self.zjump=self.values(self.l,np.arange(len(obj)))[0],self.values(self.r,np.arange(len(obj)))[0]
        edges=np.r_[self.l,self.r[-1]]
        jumps=np.r_[self.yjump[0],self.yjump[1:]-self.zjump[:-1],-self.zjump[-1]]
        if fraction>0:
            if max(abs(jumps[1:-1]))>2e-14:
                raise ValueError('Continuous reconstruction has an interior jump.')
            jumps[1:-1]=0.
        self.delta=jumps*edges/np.hypot(edges,m)
        self.jump_edges=edges
        t,w=roots_legendre(20)
        self.q=(self.l[:,None]+self.r[:,None])/2+(self.r-self.l)[:,None]*t/2
        self.w=(self.r-self.l)[:,None]*w/2
        fq,dq=self.values(self.q,np.arange(len(obj))[:,None])
        Ep=np.hypot(self.q,m)
        self.A=self.q/Ep*fq
        self.Ap=(m*m/Ep**3)*fq+self.q/Ep*dq
        self.J=float(np.sum(self.w*self.q*self.A))
        self.Jg=np.pi**2/6-quad(lambda r:r/np.expm1(r),0,cut,epsabs=1e-13)[0]

    def values(self,q,i):
        x=np.clip((q**3/3-self.U[i])/(self.V[i]-self.U[i]),0.,1.)
        d=self.z[i]-self.y[i]
        return self.y[i]+d*(x*x*(3-2*x)),d*6*x*(1-x)*q*q/(self.V[i]-self.U[i])

    def coeff(self,q,i,prime):
        f,d=self.values(q,i)
        E=np.hypot(q,self.m)
        return (self.m**2/E**3*f+q/E*d) if prime else q/E*f

    def evaluate(self,p,prime=False,ignore=None):
        # Integral over open segments plus distributional edge jumps.
        w=self.r-self.l
        near=(self.l-2*w <= p)&(p <= self.r+2*w)
        q=self.q[~near]
        Y=self.Ap if prime else self.A
        val=np.sum(self.w[~near]*Y[~near]*np.log1p(2*np.minimum(q,p)/abs(q-p)))
        for i in np.flatnonzero(near):
            l,r=self.l[i],self.r[i]
            def fun(q): return self.coeff(q,i,prime)
            if l<=p<=r:
                z=quad(lambda q:fun(q)*np.log(p+q),l,r,epsabs=2e-12)[0]
                for sg,h in ((-1,p-l),(1,r-p)):
                    if h>0:
                        z-=quad(lambda w:fun(p+sg*w),0,h,weight='alg-loga',
                                wvar=(0.,0.),epsabs=2e-12)[0]
            else:
                z=quad(lambda q:fun(q)*np.log1p(2*min(p,q)/abs(p-q)),l,r,epsabs=2e-12)[0]
            val+=z
        if prime:
            val+=2*np.sum(self.w*self.A/(p+self.q))
            use=(self.delta!=0)
            if ignore is not None:
                use[ignore]=False
            qe=self.jump_edges[use]
            val+=np.sum(self.delta[use]*np.log1p(2*np.minimum(p,qe)/abs(p-qe)))
        return float(val)

    def energy(self,p,g):
        E=np.hypot(p,self.m)
        return E+g*g*((self.Jg+self.J)/(2*np.pi**2*E)
                      -self.m**2*self.evaluate(p)/(4*np.pi**2*p*E))

    def slope(self,p,g,ignore=None):
        E=np.hypot(p,self.m)
        C=g*g*self.m**2/(4*np.pi**2*p*E)
        finite=p/E-g*g*(self.Jg+self.J)*p/(2*np.pi**2*E**3)
        return finite+C*((1/p+p/E**2)*self.evaluate(p)-self.evaluate(p,True,ignore))


class SmoothFrozen(Frozen):
    def __init__(self,d,species,m,cut,fraction):
        super().__init__(d,species,m,cut)
        self.pole=Pole(self.edges,self.f,self.m,cut,fraction)
        self.moment=self.pole.J

    @functools.lru_cache(maxsize=8192)
    def log(self,p):
        return -self.m**2*self.pole.evaluate(p)/(2*p)


def compare_support(d,grid):
    M=d['x_end']
    b=M*d['mass_ratio']
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    g=.005
    out=[]
    for cut in (.25,.5):
        rows=[]
        for fraction in (0.,.015,.05,.12):
            cls=SmoothFrozen if fraction else Frozen
            args=(fraction,) if fraction else ()
            p=cls(d,'parent',M,cut,*args)
            f=cls(d,'fermion',b,cut,*args)
            fac=solve_window(p,f,ms,Q,g*g,False)
            log=solve_window(p,f,ms,Q,g*g,True)
            rows.append(dict(transition_fraction=fraction,
                             factorized_center=fac['center'],full_center=log['center'],
                             log_center_shift=log['center']-fac['center'],
                             log_fractional_stimulus_change=log['limited_stimulus']/fac['limited_stimulus']-1,
                             full_stimulus=log['limited_stimulus']))
        baseline=rows[0]
        for row in rows:
            un=baseline['full_center']-baseline['log_center_shift']
            row['full_center_change_vs_step']=row['full_center']-baseline['full_center']
            row['factorized_center_change_vs_step']=row['factorized_center']-un
            row['full_stimulus_change_vs_step']=row['full_stimulus']/baseline['full_stimulus']-1
            row['relative_log_center_shift_change']=row['log_center_shift']/baseline['log_center_shift']-1
        out.append(dict(grid=grid,hard_cut=cut,rows=rows))
    return out


def run():
    out=[]
    support=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        d=json.load(gzip.open(path))
        support.extend(compare_support(d,grid))
        for species in ('parent','fermion'):
            edges=np.asarray(d[f'endpoint_{species}_momentum_edges'])
            f=np.asarray(d[f'endpoint_{species}_occupations'])
            m=d['x_end']*(d['mass_ratio'] if species=='fermion' else 1.)
            g=.01
            step=Pole(edges,f,m,.25,0.)
            models={x:Pole(edges,f,m,.25,x) for x in (.015,.05,.12)}
            # Probe each shell edge represented in these endpoint regions.
            for idx in np.flatnonzero((step.jump_edges>1.40)&(step.jump_edges<1.46)&(step.delta!=0.)):
                p=step.jump_edges[idx]
                C=-g*g*m*m*step.delta[idx]/(4*np.pi**2*p*np.hypot(p,m))
                finite=step.slope(p,g,ignore=idx)
                # E'=v_regular + C ln[(p+p0)/abs(p-p0)] locally.
                exponent=(-finite/C/np.log(10)) if C<0 else None
                x1,x2=p+1e-7,p+1e-11
                # Audit the log growth in the kernel's derivative itself;
                # evaluating its difference before the tiny g² is essential.
                S1=step.evaluate(x1,True)-step.evaluate(x1,True,ignore=idx)
                S2=step.evaluate(x2,True)-step.evaluate(x2,True,ignore=idx)
                predicted=step.delta[idx]*np.log((x2+p)*(x1-p)/((x1+p)*(x2-p)))
                coeff_error=abs((S2-S1)/predicted-1)
                if coeff_error>2e-12:
                    raise ValueError('Distributional shell-edge derivative failed.')
                sm=[]
                for frac,model in models.items():
                    slope=model.slope(p,g)
                    h=2e-5
                    fd=(model.energy(p+h,g)-model.energy(p-h,g))/(2*h)
                    if not slope>.99 or abs(fd-slope)>2e-10:
                        raise ValueError('Smooth slope control failed.')
                    sm.append(dict(transition_fraction=frac,finite_group_slope=slope,
                                   finite_difference_defect=float(abs(fd-slope)),
                                   real_energy_change_from_step=float(model.energy(p,g)-step.energy(p,g))))
                out.append(dict(grid=grid,species=species,edge=float(p),jump_in_occupation=float(
                                step.delta[idx]*np.hypot(p,m)/p),group_slope_log_coefficient=float(C),
                                finite_regular_group_slope=float(finite),
                                estimated_log10_inverse_distance_for_turn=exponent,
                                log_coefficient_relative_defect=float(coeff_error),smooth=sm))
    return dict(scope='Archived hard dispersion shape regularity, not a collision model.',
                derivative_coupling=.01,aperture_coupling=.005,rows=out,support=support)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
