"""Compare background stress weights on fixed archived histograms, not on an evolved relic."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import fit,scaled_energy
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack

NAMES=('count','energy','pressure','pseudo_pressure')


def weights(d,y,n):
    e=d.ds.ed
    edges=np.union1d(e,.01)
    j=np.searchsorted(e,(edges[:-1]+edges[1:])/2,side='right')-1
    a,b=edges[:-1],edges[1:]
    t,w=roots_legendre(n)
    q=a[:,None]+(b-a)[:,None]*(t+1)/2
    occ=y(j,q) if callable(y) else y[j,None]
    u=occ*(b-a)[:,None]*w*q*q/(4*np.pi**2)
    warm=np.repeat(a>=.01,n)
    return q.ravel(),u.ravel(),warm


def evaluate(q,w,warm,M):
    E=np.hypot(q,M)
    P=q*q/(3*E)
    z=np.stack((w,w*E,w*P,w*P*(q/E)**2),axis=-1)
    a=np.array([fsum(v) for v in z.T])
    b=np.array([fsum(v) for v in z[warm].T])
    return a,b


def nodal_frozen(card,base,M0):
    direct,_,_=gauss(card,base,M0,.05)
    d=ScalarGauss2(load('fine')[0],SPECS['fine'],2,4,1e-5,gauss_floor=.01)
    y=np.asarray(card['nodal_occupation'])
    n=d.ds.w*y[:d.ns]
    ends=np.r_[0,np.cumsum(d.scalar_repeats)]
    def occ(j,q):
        f=np.empty(q.shape)
        for k,i in enumerate(j):
            a,b=base.ds.ed[i:i+2]
            ss=slice(ends[i],ends[i+1])
            u=fsum(n[ss])
            if u==0:
                f[k]=0.
            elif ends[i+1]-ends[i]==1:
                f[k]=y[ends[i]]
            else:
                mu=fsum(n[ss]*scaled_energy(d.ds.p[ss],a,b,M0))/u
                r=fit(a,b,M0,u,mu)
                f[k]=np.exp(r['log_occupation_at_lower_edge']+r['tilt']*scaled_energy(q[k],a,b,M0))
        return f
    return direct,occ


def run(cards,gauss_card=None):
    l=[unpack(c,bool(i),.05) for i,c in enumerate(cards)]
    n=[c.get('warm_div',1) for c in cards]
    if n not in ([1,2],[1,2,4]):
        raise ArithmeticError('Use native then W2 and optionally W4 in that order')
    M0=l[0][0].state(.05,l[0][1])[0]
    if any(not np.array_equal(getattr(l[0][0],k),getattr(d,k)) for d,_ in l for k in ('wp','wd','qp','qd')):
        raise ArithmeticError('Charged grids changed')
    charts=[(j,weights(d,y,24),weights(d,y,48),fsum(d.ds.w*y[:d.ns])) for j,(d,y) in zip(n,l,strict=True)]
    if gauss_card is not None:
        d=l[0][0]
        expected,f=nodal_frozen(gauss_card,d,M0)
        z=[weights(d,f,j) for j in (24,48)]
        h,t=evaluate(*z[-1],M0)
        if max(abs(np.array([h[0],3*h[2],3*t[2]])/expected-1))>2e-12:
            raise ArithmeticError('Frozen positive Gauss shape does not match the accepted pressure receipt')
        charts.append(('Gauss2',*z,expected[0]))
    limiting=[]
    for level,_,sample,_ in charts:
        q,w,_=sample
        limiting.append(dict(chart=level,finite_moment_limits=[fsum(w*q**j) for j in (1,2,4)]))
    for (d,y),row in zip(l,limiting[:len(l)],strict=True):
        a,b=d.ds.ed[:-1],d.ds.ed[1:]
        exact=[fsum(y[:d.ns]*(b-a)*sum(b**(k-i-1)*a**i for i in range(k))/(2*np.pi**2*k)) for k in (4,5,7)]
        if max(abs(np.asarray(exact)/row['finite_moment_limits']-1))>2e-14:
            raise ArithmeticError('Polynomial extreme-mass moment test failed')
    out=[]
    for M in (0.,1e-6,1e-5,float(M0),.01,.1,1.,3.,10.,100.):
        A=[]
        for level,low,hi,count in charts:
            a,b=evaluate(*hi,M)
            lo=evaluate(*low,M)
            err=max(np.max(abs(x/z-1)) for x,z in zip(lo,(a,b),strict=True))
            if err>2e-12 or not 0<a[3]<=a[2]<=a[1]/3 or abs(a[0]/count-1)>3e-14:
                raise ArithmeticError('Frozen stress moments did not close')
            span=1e-4
            lag=[evaluate(*hi,M*np.exp(t))[0]*np.exp(np.array([-3,-4,-4,-4])*t) for t in (-span,span)]
            slope=(lag[1]-lag[0])/(2*span)
            ref=np.array([-3*a[0],-3*(a[1]+a[2]),-5*a[2]+a[3]])
            ident=max(abs(slope[:3]/ref-1))
            if ident>1e-7:
                raise ArithmeticError('Independent rescaling derivative check failed')
            p=a[2]/a[1]
            A.append(dict(warm_chart=level,totals=a.tolist(),tail_fraction=(b/a).tolist(),p_over_rho=p,
                          pseudo_pressure_over_p=a[3]/a[2],
                          adiabatic_derivative=p*(5-a[3]/a[2])/(3*(1+p)),
                          max_paired_quadrature_defect=err,max_dilation_derivative_defect=ident))
        out.append(dict(test_mass=M,rows=A,changes_relative_to_native=[
            (np.asarray(r['totals'])/np.asarray(A[0]['totals'])-1).tolist() for r in A[1:]]))
    raw=np.asarray([v['finite_moment_limits'] for v in limiting])
    if np.max(abs(np.asarray([r['totals'][1] for r in out[0]['rows']])/raw[:,0]-1))>2e-14:
        raise ArithmeticError('The exact massless density limit failed')
    return dict(fixed_evolved_x=.05,tail_cut=.01,reference_fitted_mass=M0,columns=NAMES,measurements=out,
                frozen_limit_columns=('massless_energy','heavy_limit_3MP','heavy_limit_3M3PseudoPressure'),
                limits=limiting,relative_limits_to_native=(raw[1:]/raw[0]-1).tolist(),
                scope='Vary a formal M=a*m/T at fixed archived f for stress-weight diagnostics only. No clock/scale-factor assignment or physical late-time transfer.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path,nargs='+')
    p.add_argument('--gauss',type=Path)
    a=p.parse_args()
    print(json.dumps(run([read(f) for f in a.state],read(a.gauss) if a.gauss else None),indent=2,allow_nan=False))
