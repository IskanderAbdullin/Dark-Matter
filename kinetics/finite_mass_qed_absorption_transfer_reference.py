"""Photon-energy-first absorption using both timelike and spacelike neutral transfer."""
import gzip,json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q
from early_pole_hard_yukawa_photoproduction_reference import invariant_form
from finite_mass_qed_signed_support_bridge_reference import BridgeSource


@lru_cache(None)
def gl(n):return roots_legendre(n)


def transfer(M,m,s,k,c):
    Es=np.hypot(s,Q)
    J=k*(Es-Q*c)
    U=s*s-2*J
    E=Es-k
    bz=Q*c-k
    st=Q*np.sqrt((1-c)*(1+c))
    V=np.hypot(bz,st)
    D=M*M+U-m*m
    H=np.sqrt(D*D-4*M*M*U)
    upper_expression=E*D+V*H
    lo=(D*D+4*M*M*V*V)/(2*upper_expression)
    high=upper_expression/(2*U) if U>0 else np.inf
    return Es,J,U,E,V,D,lo,high,bz/V,st/V


def average(M,m,k,J,U,E,V,D,c,s,ep):
    ep=np.asarray(ep)
    eb=ep-E
    pz=(E*ep-D/2)/V
    pbz=pz-V
    pp2=np.maximum(0.,(ep-M)*(ep+M)-pz*pz)
    bmag=abs(pbz)
    stable=(m*m+pp2)/(eb+bmag)+bmag*s*s/(1+abs(c))
    C=k*np.where(pbz*c>0.,stable,eb-pbz*c)
    A=C+J
    yt=k*np.sqrt(pp2)*s
    Dp=np.sqrt((A-yt)*(A+yt))
    Db=np.sqrt((C-yt)*(C+yt))
    d=(M*M+m*m-U)/2
    q0=4*(J*J-2*J*d+2*d*M*m+2*d*d)
    val=q0*(A+C)/((Dp+Db)*Dp*Db)-4*(-J+M*m+d)*(M*M*A/Dp**3+m*m*C/Db**3)
    return val,A,C,yt


def edge_U(M,m,e,E):
    p=np.sqrt((E-M)*(E+M))
    d=2*e*E-M*M+m*m
    center=d-2*p*p
    eb=E-e
    delta=2*p*np.sqrt((eb-m)*(eb+m))
    l=center+np.copysign(delta,center)
    prod=(2*e*M*M/(E+p)-M*M+m*m)*(2*e*(E+p)-M*M+m*m)
    h=np.divide(prod,l,out=center.copy(),where=l!=0.)
    return np.r_[l,h]


def angle_pilot(M,m,s):
    # Time-routed transfer changes character near k=s²/(2Es).
    errs=[]; endpoints=[]
    x,w=gl(240)
    for k in (4e-5,2e-4,.0007,.0013,.0025):
        for c in (-.875,-.31,.04,.77,.99):
            Es,J,U,E,V,D,l,h,co,si=transfer(M,m,s,k,c)
            l=max(l,M,E+m)
            h=min(h,22.)
            if h<=l:continue
            for t in (.02,.32,.83):
                z=np.array([l+t*(h-l)])
                v,A,C,b=average(M,m,k,J,U,E,V,D,co,si,z)
                y=b[0]*np.cos(np.pi*(x+1))
                exact=w@invariant_form(A[0]-y,C[0]-y,(M*M+m*m-U)/2,M,m)/2
                errs.append(float(abs(v[0]/exact-1.)))
            if l>max(M,E+m)+1e-11:
                p=np.sqrt((l-M)*(l+M))
                endpoints.append(abs(abs((E*l-D/2)/(V*p))-1.))
    if max(errs)>3e-8 or max(endpoints)>3e-10:
        raise ArithmeticError((max(errs),max(endpoints)))
    return dict(direct_azimuth_quadrature_max_relative_defect=max(errs),
                analytic_lower_parent_endpoint_max_cosine_defect=float(max(endpoints)))


def point(d,k,nc,ne):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es=np.hypot(s,Q)
    Pe=np.hypot(M,d['endpoint_parent_momentum_edges'])
    Be=np.hypot(m,d['endpoint_fermion_momentum_edges'])
    fp,fb=[np.r_[0.,d[f'endpoint_{z}_occupations'],0.] for z in ('parent','fermion')]
    e=Es-k
    e_edges=np.unique(np.r_[Pe,Be+e])
    ok=(e_edges>M)&(e_edges>e+m)
    Us=edge_U(M,m,e,e_edges[ok])
    Cs=(Es+(Us-s*s)/(2*k))/Q
    split=np.unique(np.r_[-1.,-.999,-.985,-.91,-.65,-.25,.25,.65,.91,.985,.999,1.,
                         Cs[(Cs>-1.)&(Cs<1.)]])
    cz,wz=gl(nc);z,w=gl(ne)
    total=np.zeros(2)
    for a,b in zip(split[:-1],split[1:],strict=True):
        for x,v in zip(cz,wz,strict=True):
            C=(a+b+(b-a)*x)/2
            _,J,U,e,V,D,l,h,co,si=transfer(M,m,s,k,C)
            bounds=(l,h)
            l=max(l,M,e+m)
            h=min(h,e_edges[-1])
            if h<=l:continue
            guides=[]
            ladd=np.array([.006,.032,.15,.7,3.,14.,65.,300.,1500.,7000.])
            for ep,sg in zip(bounds,(1.,-1.),strict=True):
                if not np.isfinite(ep) or not (l-.001*(h-l)<=ep<=h+.001*(h-l)):continue
                pz=(e*ep-D/2)/V
                deriv=abs(2*(ep-e*pz/V))
                for mass in (M,m):
                    dx=mass*mass/max(deriv,1e-15)
                    guides.extend(ep+sg*dx*ladd)
            cU=e*e-(Q*C-k)**2
            # The collinear case excites small charged denominators. Split at
            # the roots of the same charged energies in an axis along k.
            for mass,A,shift in ((M,D/2,0.),(m,D/2-U,e)):
                disc=A*A-mass*mass*cU
                if disc>0.:
                    big=e*A+np.copysign(abs(Q*C-k)*np.sqrt(disc),e*A)
                    for q in (big/cU,(A*A+mass*mass*(Q*C-k)**2)/big):
                        zz=q+shift
                        if l-.03*(h-l)<zz<h+.03*(h-l):
                            pz=(e*zz-D/2)/V
                            grad=abs(2*(zz-e*pz/V))
                            scale=(mass*mass+2*mass*abs(pz*si))/max(grad,1e-15)
                            guides.extend(zz+scale*np.r_[-ladd,0.,ladd])
            cuts=np.unique(np.clip(np.r_[l,h,e_edges,guides],l,h))
            l,h=cuts[:-1],cuts[1:]
            Y=(l+h)[:,None]/2+(h-l)[:,None]*z/2
            A,*_=average(M,m,k,J,U,e,V,D,co,si,Y)
            if A.min()<-1e-8:raise ArithmeticError('Negative absorptive matrix element')
            med=(l+h)/2
            F=fp[np.searchsorted(Pe,med)]
            B=fb[np.searchsorted(Be,med-e)]
            W=(h-l)/2*(A@w)
            total+=v*(b-a)/2*k/(128*np.pi**3*Es*V)*np.array([W@(F*(1-B)),W@((1-F)*B)])
    return total


def run(grids):
    root=Path(__file__).parent/'results'
    ms=((19,44),(37,58))
    out={}
    for g in grids:
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        zz=[]
        for k in (.00016,.0002,.0004,.000435,.00045,.0008,.00128,.00256):
            a,b=(point(d,k,*mesh) for mesh in ms)
            zz.append(dict(k=k,gain_loss_before_Planck=b.tolist(),two_mesh_relative_drift=float(max(abs(a/b-1.)))))
        out[g]=zz
    return dict(meshes=[list(x) for x in ms],data=out)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('third','fine')),indent=2,allow_nan=False))
