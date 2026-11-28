"""Massive background and derivative moment consistency on candidate CLASS grids."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre

T0=8.617333262e-5*2.72548
A_ARCH=5.326008998418959e-10


def spline_tail(path):
    q,f=np.loadtxt(path).T
    assert np.all(f>0) and np.all(np.diff(q)>0)
    F=CubicSpline(q,f)
    decay=(f[-2]-f[-1])/(q[-1]-q[-2])/f[-1]
    assert decay>0
    def value(p):
        y=np.clip(p,q[0],q[-1])
        t=f[-1]*np.exp(-decay*np.maximum(p-q[-1],0.))
        z=np.where(p>q[-1],t,F(y))
        dz=np.where(p>q[-1],-decay*t,np.where(p<q[0],0.,F(y,1)))
        assert np.all(z>=0)
        return z,dz
    return q,value


def moments(q,w,f,mT):
    y,d=f(q)
    E=np.hypot(q,mT)
    R=np.dot(w,q*q*E*y)
    P=np.dot(w,q**4*y/(3*E))
    J=np.dot(w,q**6*y/(3*E**3))
    return np.array([R,P]),np.array([np.dot(w,q**3*E*d),np.dot(w,q**5*d/(3*E))]),J


def boundary(x,f,mT):
    if x is None:
        return np.zeros(2)
    E=np.hypot(x,mT)
    h=float(f(np.asarray(x))[0])
    return np.array([x**3*E,x**5/(3*E)])*h


def benchmark():
    f=lambda q:(np.exp(-q),-np.exp(-q))
    x,w=roots_legendre(20)
    v=np.r_[0.,np.geomspace(1e-6,.9,42),np.arange(1.,11.)]
    q=((v[:-1,None]+v[1:,None])+(v[1:,None]-v[:-1,None])*x)/2
    w=(w*(v[1:,None]-v[:-1,None])/2).ravel()
    ds=[]
    for y in (.01,1.,100.):
        a,d,j=moments(q.ravel(),w,f,y)
        E=boundary(10.,f,y)-np.array([3*(a[0]+a[1]),5*a[1]-j])
        ds.append(np.max(abs(d/E-1.)))
    assert max(ds)<1e-12
    return max(ds)


def native(n,qmax):
    if qmax is not None:
        h=qmax/n
        q=h*np.arange(1,n+1)
        w=np.full(n,h)
        w[-1]*=.5
    else:
        t=np.arange(1,n+1)/(n+1)
        q=1/t-1
        w=1/(n+1)/t**2
    return q,w


def continuum(knots,qmax):
    upper=(2*knots[-1]+20.) if qmax is None else qmax
    e=np.unique(np.clip(np.r_[0.,knots,knots[-1]+np.geomspace(1e-3,upper,80),upper],0.,upper))
    x,w=roots_legendre(12)
    h=(e[1:]-e[:-1])[:,None]/2
    q=(e[1:]+e[:-1])[:,None]/2+h*x
    return q.ravel(),(h*w).ravel()


def case(path,m,n,qmax,label):
    ks,f=spline_tail(path)
    p,w=native(n,qmax)
    q0,w0=continuum(ks,qmax)
    tau=m/T0
    a=np.sort(np.r_[[1e-13,1e-11,A_ARCH,1e-7,1e-4,1e-2,1.],
                    np.array([.1,1.,4.,16.,64.])/tau])
    rows=[]
    for scale in a:
        x=scale*tau
        u,d,j=moments(p,w,f,x)
        v,e,jr=moments(q0,w0,f,x)
        b=boundary(qmax,f,x)
        dt=b-np.array([3*(u[0]+u[1]),5*u[1]-j])
        rt=b-np.array([3*(v[0]+v[1]),5*v[1]-jr])
        err=e/rt-1.
        assert max(abs(err))<1e-8
        rows.append(dict(scale_factor=float(scale),am_over_Tgamma0=float(x),
                energy_and_pressure_quadrature_error=(u/v-1.).tolist(),
                energy_and_pressure_dilation_residual=(d/dt-1.).tolist(),
                continuum_crosscheck=err.tolist()))
    e=np.asarray([k["energy_and_pressure_quadrature_error"] for k in rows])
    d=np.asarray([k["energy_and_pressure_dilation_residual"] for k in rows])
    return dict(label=label,mass_eV=m,momentum_nodes=n,qmax=qmax,rows=rows,
                max_abs_energy_and_pressure_quadrature_error=np.max(abs(e),axis=0).tolist(),
                max_abs_energy_and_pressure_dilation_residual=np.max(abs(d),axis=0).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("fine",type=Path);p.add_argument("stable_balanced",type=Path);p.add_argument("scalar_balanced",type=Path)
    a=p.parse_args()
    s=a.fine/"stable.dat";v=a.fine/"scalar.dat"
    cfg=[(s,2e7,n,20.,f"stable_raw{n}") for n in (64,128,256)]
    cfg.extend([(s,2e7,128,21.370224719101124,"stable_raw128_qshift"),
                (a.stable_balanced,2e7,256,20.,"stable_smooth256"),
                (s,2e7,512,None,"stable_raw512_inverse_q_mapping"),
                (v,86.,160,63.8376753507014,"scalar_raw160_qshift"),
                (v,86.,384,96.,"scalar_raw384"),
                (a.scalar_balanced,86.,128,96.,"scalar_smooth128"),
                (a.scalar_balanced,86.,256,96.,"scalar_smooth256")])
    print(json.dumps(dict(synthetic_nonzero_boundary_max_defect=benchmark(),
              cases=[case(*row) for row in cfg],
              scope="External numerical audit with the tabulated spline's analytic derivative. Native nodes and weights, finite boundary term included. Not native perturbation or a physical transfer estimate."),
              indent=2,allow_nan=False))
