"""Exact hard 2↔2 subcell triangle with first-order reciprocal dispersions."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre,expit

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon
from early_pole_gauge_projection_momentum_reference import curve


def inv(E,m,A):
    q=(E+np.sqrt(np.maximum(E*E-4*A,0.)))/2
    return np.sqrt(np.maximum((q-m)*(q+m),0.))


def pair(d,g,label,case,order):
    q1,q2,ct=case
    x=d['x_end']
    mass=x if label=='parent' else d['mass_ratio']*x
    J=[]
    for k,m in [('parent',x),('fermion',d['mass_ratio']*x)]:
        J.append(cell_J(d[f'endpoint_{k}_momentum_edges'],
                        d[f'endpoint_{k}_occupations'],m,.5))
    jg=photon(1.,-.2,.5)
    A=g*g/np.pi**2*sum(J)
    B=.5*g*g/np.pi**2*(jg+J[label=='fermion'])
    P=np.sqrt(q1*q1+q2*q2+2*q1*q2*ct)
    E=curve(q1,0.,A)+curve(q2,0.,A)
    pp=inv(E-curve(.5,mass,B),mass,B)
    if pp<=.5:return None
    f=lambda q:inv(E-curve(q,mass,B),mass,B)
    axes=[lambda q:q-f(q)-P,lambda q:f(q)-q-P,lambda q:P-q-f(q)]
    e=np.asarray(d[f'endpoint_{label}_momentum_edges'])
    nodes=np.linspace(.5,pp,76)
    sep=np.r_[.5,pp]
    for h in axes:
        v=h(nodes)
        for i in np.flatnonzero(v[:-1]*v[1:]<0):
            sep=np.r_[sep,brentq(h,nodes[i],nodes[i+1],xtol=3e-14)]
    # Sample the physical aperture before splitting for any shell changes.
    ok=[(l,r) for l,r in zip(np.sort(sep)[:-1],np.sort(sep)[1:])
        if max(h((l+r)/2) for h in axes)<=0.]
    if not ok:return None
    if g==0.:
        s=E*E-P*P
        lo=.5*(E-P*np.sqrt(1-4*mass*mass/s))
        hi=.5*(E+P*np.sqrt(1-4*mass*mass/s))
        ref=[max(.5,np.sqrt(lo*lo-mass*mass)),
             min(pp,np.sqrt(hi*hi-mass*mass))]
        assert max(abs(np.array([ok[0][0],ok[-1][1]])-ref))<3e-12
    kinks=np.r_[e,inv(E-curve(e,mass,B),mass,B)]
    z,w=roots_legendre(order)
    p,k,ww=[],[],[]
    for lo,hi in ok:
        tt=np.unique(np.r_[lo,kinks[(kinks>lo)&(kinks<hi)],hi])
        th=np.arcsin(np.clip(2*(tt-lo)/(hi-lo)-1.,-1.,1.))
        dd=np.diff(th)[:,None]/2
        t=(th[1:]+th[:-1])[:,None]/2+dd*z
        u=(lo+(hi-lo)*(.5+.5*np.sin(t))).ravel()
        p.extend(u);k.extend(f(u));ww.extend((dd*w*np.cos(t)*(hi-lo)/2).ravel())
    p,k,w=np.array(p),np.array(k),np.array(ww)
    eP,eB=curve(p,mass,B),curve(k,mass,B)
    mu1=(P*P+q1*q1-q2*q2)/(2*P*q1)
    mu=(P*P+p*p-k*k)/(2*P*p)
    assert max(abs(mu))<=1.
    pho1=np.r_[q1*np.sqrt(1-mu1*mu1),q1*mu1]
    pho2=np.r_[0.,P]-pho1
    a=np.column_stack([p*np.sqrt(1-mu*mu),p*mu])
    b=np.r_[0.,P]-a
    # A positive geometry weight, not a gauge spin-dependent matrix element.
    e0=np.hypot(k,mass)
    vel=(1-B/(e0*e0))*k/e0
    w=w*p*k/(P*eP*eB*vel)
    i,j=(np.searchsorted(e,v)-1 for v in (p,k))
    F=np.asarray(d[f'endpoint_{label}_occupations'])
    g1,g2=1/np.expm1(np.array([curve(q1,0.,A),curve(q2,0.,A)])+.2)
    forward=g1*g2*(1-F[i])*(1-F[j])
    inverse=(1+g1)*(1+g2)*F[i]*F[j]
    r=w*(forward-inverse)
    z=np.zeros((len(F),4))
    np.add.at(z,i,np.column_stack([r,r*eP,r[:,None]*a]))
    np.add.at(z,j,np.column_stack([r,r*eB,r[:,None]*b]))
    R=z.sum(axis=0)-r.sum()*np.r_[2.,E,pho1+pho2]
    assert np.max(abs(R))/max(sum(abs(r)),1e-40)<8e-15
    # Separate detailed-balance audit at common T and common species chemical
    # potential for the weighted count N_gamma+sum N_fermion.
    fa,fb=expit(-eP-.2),expit(-eB-.2)
    eq=g1*g2*(1-fa)*(1-fb)/((1+g1)*(1+g2)*fa*fb)
    assert max(abs(eq-1))<5e-14
    return dict(g=g,species=label,incoming=case,momentum_bounds=ok,
                frozen_signed_geometry_weighted_flux=r.sum(),phase_weight=w.sum(),
                total_count_energy_vector_defect=(R/max(sum(abs(r)),1e-40)).tolist(),
                thermal_relative_defect=max(abs(eq-1)))


def run():
    d=load('results/early-pole-spectrum-third-m06-2026-11-26.json.gz')
    rows=[]
    for case in ((.72,1.25,-.25),(.55,1.93,.65),(1.08,2.32,.92),(1.4,1.3,-.92)):
        for g in (0.,.005,.05):
            for l in ('parent','fermion'):
                a,b=(pair(d,g,l,case,n) for n in (12,28))
                if a is None:continue
                err=abs(a['phase_weight']/b['phase_weight']-1)
                assert err<2e-10
                rows.append(dict(phase_refinement=err,**b))
    return dict(scope='Fixed 2↔2 continuous hard event geometry and reversible occupancies. '
                      'No coupling coefficient, scattering spin amplitude or evolving photon bath.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
