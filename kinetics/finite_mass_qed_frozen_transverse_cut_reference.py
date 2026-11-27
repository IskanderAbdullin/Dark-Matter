"""Frozen isotropic transverse Landau cut and pole spectral sum rule."""
import gzip,json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF,Q,Species,case,transverse


def response(charges,u):
    real=0.
    signed=0.
    positive=0.
    for c in charges:
        # The response has been integrated by parts in velocity for the
        # imaginary side, giving a positive integral even for histograms.
        q=u/np.sqrt((1-u)*(1+u))
        p0=c.m*q
        lower=np.maximum(c.e[:-1],p0)
        each=np.maximum((c.e[1:]-lower)*(c.e[1:]+lower),0.)
        positive+=np.pi*u*(1-u*u)*(c.f@each)/4
        mask=(c.v>u)
        z=u/c.v[mask]
        assert np.all(z<1.)
        h=.5*(z*z+z*(1-z*z)*np.arctanh(z))
        real+=((c.W[mask]*c.v[mask]**2)@h)
        signed+=np.pi/4*((c.W[mask]*c.v[mask]**2)@(z*(1-z*z)))
        tail=(~mask)&(c.v>0)
        x=c.v[tail]/u
        # The real value at the integrable cusp z=1 is exactly 1/2.
        r=np.full_like(x,.5)
        p=x<1.
        r[p]=transverse(x[p])[0]
        real+=((c.W[tail]*c.v[tail]**2)@r)
    return PREF*real,PREF*positive,PREF*signed


def validate(d):
    M=d['x_end'];m=M*d['mass_ratio']
    es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    ch=[Species(x,d[f'endpoint_{name}_momentum_edges'],d[f'endpoint_{name}_occupations'])
        for x,name in ((M,'parent'),(m,'fermion'))]
    poles=case(d)['poles']
    angles=np.unique(np.concatenate([c.v for c in ch]))
    angles=angles[(angles>0.)&(angles<1.)]
    test=[response(ch,x) for x in (.00001,.0001,.003,.031,.27,.64,.94,.994,.9999,.999999)]
    scale=sum(c.D for c in ch)
    ident=max(abs(a-b) for _,a,b in test)/scale
    if ident>6e-14 or min(a for _,a,b in test)<0.:
        raise ArithmeticError('The retarded Landau cut sign or infrared primitive failed')
    tracks=[]
    for obj in poles:
        k=obj['k']
        if k==0.:continue
        Z=obj['physical_transverse_w2_pole_residue']
        def integrand(u):
            r,i,_=response(ch,u)
            A=-k*k*(1-u*u)-r
            rho=2*i/(A*A+i*i)
            return k*k*u*rho/np.pi
        def interval(u0,u1,tol):
            if u1<=u0:return 0.
            points=np.unique(np.r_[angles,np.geomspace(1e-12,.1,18),.95,.99,.999,.9999]) 
            points=points[(points>u0)&(points<u1)]
            return quad(integrand,u0,u1,epsabs=tol,epsrel=2e-9,points=points,limit=len(points)+160)[0]
        I,Ih=(interval(0.,1.,tol) for tol in (2e-10,2e-12))
        W=Z+Ih
        if abs(W-1.)>8e-9:
            raise ArithmeticError(('Transverse cut plus pole sum rule',k,Z,Ih,W))
        cold=min(es/k,1.)
        frac=Ih if cold==1. else interval(0.,cold,2e-12)
        tracks.append(dict(k=k,pole_weight=Z,landau_cut_weight=Ih,
                           canonical_weight_defect=W-1.,mesh_cut_abs_change=abs(I-Ih),
                           cut_weight_below_scalar_energy=frac))
    support=[]
    for k,u,c in ((.0008,.25,0.),(.00128,.5,-.5),(.00256,.25,.5)):
        omega=k*u
        E=es-omega
        V=np.hypot(k-Q*c,Q*np.sqrt((1-c)*(1+c)))
        U=(E-V)*(E+V)
        D=M*M+U-m*m
        H=np.sqrt(D*D-4*M*M*U)
        Emin=(D*D+4*M*M*V*V)/(2*(E*D+V*H))
        Emin=max(Emin,M,m+E)
        # Locate a sample parent energy inside the actual stored occupied
        # intervals for both charged spectra, not merely below the caps.
        ep=Emin*1.01
        pp,bp=np.sqrt((ep-M)*(ep+M)),np.sqrt((ep-E-m)*(ep-E+m))
        f=[sp.f[np.searchsorted(sp.e,w)-1] for sp,w in zip(ch,(pp,bp),strict=True)]
        R,I,_=response(ch,u)
        if U>=0 or omega>=es or I<=0 or not all(x>0 for x in f):
            raise ArithmeticError(('Manufactured support witness failed',k,u,c,U,Emin,f))
        support.append(dict(k=k,omega=omega,cosine=c,minimum_parent_energy=Emin,
                            nonzero_parent_and_daughter_fixture_energy=ep,occupations=f,
                            minus_retarded_transverse_Im=I))
    return dict(positive_radial_identity_max_defect_in_Debye_units=ident,rows=tracks,
                below_scalar_energy_positive_cut_with_charged_support=support)


def run():
    root=Path(__file__).parent/'results'
    out=[]
    for g in ('fine','third'):
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        out.append(dict(grid=g,**validate(d)))
    return dict(scope='Collisionless isotropic transverse spectral density, no longitudinal conversion or nonthermal photon statistical law.',data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
