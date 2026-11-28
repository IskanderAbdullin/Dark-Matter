"""Frozen longitudinal Coulomb response: positive cut and static sum rule."""

import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from finite_mass_qed_frozen_longitudinal_pole_reference import L
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF, Species
from finite_mass_qed_spacelike_density_statistics_reference import exact, soft


class DensityResponse:
    def __init__(self, charges):
        self.charges=charges
        pv=np.concatenate([c.v for c in charges])
        pw=np.concatenate([c.W for c in charges])
        self.v,ix=np.unique(pv,return_inverse=True)
        weights=np.bincount(ix,weights=pw)
        keep=(weights!=0)&(self.v>0)
        self.W,self.v=weights[keep],self.v[keep]
        self.debye=PREF*sum(self.W)

    def retarded(self,u):
        v,W=self.v,self.W
        near=(v>u)
        low=(v<.2*u)
        reg=~(near|low)
        if np.any(v==u):
            return np.nan,np.nan
        M=np.empty_like(v)
        M[low]=L(v[low]/u)[0]
        M[near]=(u/v[near])*np.arctanh(u/v[near])-1
        M[reg]=(u/v[reg])*np.arctanh(v[reg]/u)-1
        R=PREF*(W @ M)
        I=PREF*np.pi*u/2*(W[near] @ (1/v[near]))
        return R,I

    def direct_angular(self,u):
        out=0.
        for v,w in zip(self.v,self.W,strict=True):
            z=u/v
            if z<1.:
                h=-quad(lambda t:t/2,-1,1,weight='cauchy',wvar=z,
                        epsabs=2e-13,epsrel=2e-13)[0]
            else:
                h=quad(lambda t:t/(2*(z-t)),-1,1,epsabs=2e-13,epsrel=5e-12)[0]
            out+=w*h
        return PREF*out

    def pole(self,k):
        u=brentq(lambda x:self.retarded(x)[0]-k*k,1.,max(2.,.04/k),xtol=6e-15)
        pi_prime= -PREF*(self.W*self.v @ L(self.v/u)[1])/(k*u*u)
        return k*u,-1/pi_prime

    def samples(self,u):
        U=np.asarray(u)[:,None]
        V=self.v[None,:]
        dif=abs(U-V)
        with np.errstate(divide='ignore',invalid='ignore'):
            Z=V/U
            C=.5*(U/V)*np.log1p(2*np.minimum(U,V)/dif)-1
            low=Z<.2
            C[low]=L(Z[low])[0]
            C[dif==0.]=0.
        R=PREF*(C@self.W)
        I=PREF*np.pi*u/2*((U<V)@(self.W/self.v))
        R[np.any(dif==0.,axis=-1)]=np.inf
        if np.min(I)<-1e-18:
            raise ArithmeticError('Imaginary density response violates isotropic positivity')
        return R,np.maximum(I,0.)


def quadrature(model,endpoints,cut,order):
    knots=np.unique(np.clip(np.r_[endpoints,0.,cut],0.,cut))
    x,w=roots_legendre(order)
    theta=np.pi*(x+1)/4
    left,right=knots[:-1,None],knots[1:,None]
    nodes=(left+(right-left)*np.sin(theta)**2).ravel()
    weight=((right-left)*(np.pi/4)*w*np.sin(2*theta)).ravel()
    R,I=model.samples(nodes)
    return R,I,weight/(np.pi*nodes),nodes


def investigate(g,d):
    M=d['x_end'];m=M*d['mass_ratio']
    cs=[Species(v,d[f'endpoint_{n}_momentum_edges'],d[f'endpoint_{n}_occupations'])
        for v,n in ((M,'parent'),(m,'fermion'))]
    eforms=[(c.m,c.E,c.f) for c in cs]
    model=DensityResponse(cs)
    u0=(.0047,.057,.25,.67,.9,.99,.9996)
    checks=[]
    for u in u0:
        real,imag=model.retarded(u)
        z=np.sum([soft(*x,u) for x in eforms],axis=0)
        indep=PREF*np.pi*u/2*z[1]
        checks.append(abs(imag/indep-1))
    realcheck=max(abs(model.retarded(u)[0]-model.direct_angular(u))/model.debye for u in (.003,.17,.67,.91))
    if realcheck>2e-10 or max(checks)>2e-12:
        raise ArithmeticError((g,realcheck,checks))
    # All of these modes have a timelike pole and no real pole beyond the
    # archived Landau support before u=1.
    K=(.0002,.0008,.00128,.00256,.004)
    # Angular jump locations are integrable logarithmic features; split them explicitly.
    pts=np.unique(np.r_[1e-9,.01,.1,.9,.99,.999,.9999,model.v])
    pts=pts[(pts>0)&(pts<1)]
    results=[]
    qlo=quadrature(model,pts,1.,36)
    qhi=quadrature(model,pts,1.,128)
    for k in K:
        pole,res=model.pole(k)
        def folded(q,h=0.):
            R,I,W,u=q
            if h:
                W=W*u*u/(u*u+h*h)
            return float((W @ (2*I/((k*k-R)**2+I*I)))*k*k)
        cut=folded(qhi)
        independent=folded(qlo)
        polepart=2*k*k*res/pole
        target=model.debye/(k*k+model.debye)
        cold=min(d['checkpoints'][-1]['scalar_pole_over_T']/k,1.)
        low=cut if cold==1. else folded(quadrature(model,pts,cold,96))
        euclidean=[]
        for h in (.03,.3,3.):
            z=model.v/h
            Pi=PREF*(model.W@(np.arctan(z)/z-1))
            direct=-Pi/(k*k-Pi)
            spectral=folded(qhi,h)+polepart*pole*pole/(pole*pole+(k*h)**2)
            euclidean.append(abs(spectral-direct))
        err=max(abs(cut-independent),abs(cut+polepart-target),max(euclidean))
        if min(cut,polepart,target)<=0 or not 0<low<=cut+5e-11 or err>5e-9:
            raise ArithmeticError((g,k,cut,polepart,target,err,low))
        results.append(dict(k=k,pole=pole,residue_in_omega=res,cut_normalized=cut,pole_normalized=polepart,
                            static_debye_target=target,static_sum_defect=cut+polepart-target,
                            independent_order_cut_defect=abs(cut-independent),imaginary_axis_max_defect=max(euclidean),
                            cut_below_cold_scalar_energy=low))
    finite=[]
    for k in (.0002,.0008,.00256,.01):
        for u in u0:
            w=k*u
            exact_raw=np.sum([exact(*x,k,w) for x in eforms],axis=0)[3]
            finite.append(abs((PREF*np.pi*w*exact_raw/(16*k))/model.retarded(u)[1]-1))
    return dict(grid=g,real_cauchy_control_defect_in_debye_units=realcheck,
                soft_density_retarded_ibp_max_relative=max(checks),finite_free_vs_soft_imaginary_max_relative=max(finite),
                data=results)


if __name__=='__main__':
    p=Path(__file__).parent/'results'
    data=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as src:
            data.append(investigate(g,json.load(src)))
    print(json.dumps(dict(scope='Frozen one-loop Coulomb pole plus positive soft cut, no occupation of collective modes.',data=data),
                     indent=2,allow_nan=False))
