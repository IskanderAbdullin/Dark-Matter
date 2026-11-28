"""A bounded spacelike current-noise probe with no charged intermediate roots."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I,dot,slash,norm_spinors
from finite_mass_qed_frozen_transverse_cut_reference import response
from finite_mass_qed_longitudinal_static_spectral_reference import DensityResponse
from finite_mass_qed_spacelike_density_statistics_reference import soft
from finite_mass_qed_transverse_kinetic_noise_reference import archive
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup
from finite_mass_qed_timelike_pole_batched_absorption_reference import load


def spectral(ch,k,w,name):
    u=w/k
    if name=='T':
        R,S,_=response(ch,u)
        a=w*w-k*k-R
        r=np.array([archive(c,u,48) for c in ch])
    else:
        R,S=DensityResponse(ch).retarded(u)
        a=k*k-R
        r=np.array([soft(c.m,c.E,c.f,u) for c in ch])
    noise=r[:,0].sum()/r[:,1].sum()
    rho=2*S/(a*a+S*S)
    large=rho*(noise/w+.5)
    small=rho*(noise/w-.5)
    if not 0.<small<large or abs((large-small)/rho-1.)>3e-12:
        raise ArithmeticError('A positive continuum spectral split is not reconstructed')
    return small,large,rho,noise


def envelopes(M,m,t,V,az,ax,k,w,s0,qK,lo,hi):
    scalar_prod=(M*M+(t-V)*(t+V)-m*m)/2
    def x(E):return (t*E-scalar_prod)/V
    def r(E):return np.sqrt((E-M)*(E+M)-x(E)**2)
    B=2*k*abs(ax)
    dc=2*(w-k*az*t/V)
    def diff(E,h):return dc+h*B*(E-x(E)*t/V)/r(E)
    def zz(E,h):
        return s0+2*(w*E-k*az*x(E))+h*B*r(E)
    Q=[]
    for h in (-1.,1.):
        vals=[zz(lo,h),zz(hi,h)]
        if diff(lo,h)*diff(hi,h)<0.:
            E=brentq(lambda v:diff(v,h),lo,hi)
            vals.append(zz(E,h))
        Q.extend(vals)
    A,B=min(Q),max(Q)
    res=[(A,B),(2*qK-B,2*qK-A)]
    if any(a*b<=0 for a,b in res):raise ArithmeticError(('A patch has an on-shell line',res))
    return res,min(abs(x) for z in res for x in z)


def panel(d,ch,k,w,c,name,limits,n,az):
    M,m,q,K,T,_=setup(d,k,w,c)
    V=np.linalg.norm(T[1:]);t=T[0]
    axis=T[1:]/V
    cross=np.array([axis[2],0.,-axis[0]])
    lo,hi=limits
    s0=dot(K,K)
    signs,gap=envelopes(M,m,t,V,axis[2],axis[0],k,w,s0,dot(q,K),lo,hi)
    ed=np.unique(np.clip(np.r_[lo,hi,ch[0].E,ch[1].E+t],lo,hi))
    a,b=ed[:-1],ed[1:];mid=(a+b)/2
    f,g=[np.r_[0.,c.f,0.][np.searchsorted(c.E,j)] for c,j in zip(ch,(mid,mid-t),strict=True)]
    x,h=roots_legendre(n)
    E=(mid[:,None]+(b-a)[:,None]*x/2).ravel()
    W=((b-a)[:,None]*h/2).ravel()
    X=(t*E-(M*M+dot(T,T)-m*m)/2)/V
    C=np.sqrt((E-M)*(E+M)-X*X)
    phi=(np.arange(az)+.5)*2*np.pi/az
    tree=np.cos(phi)[:,None]*cross+np.sin(phi)[:,None]*np.array([0.,1.,0.])
    xyz=X[:,None,None]*axis+C[:,None,None]*tree
    p=np.concatenate([np.broadcast_to(E[:,None,None],(*xyz.shape[:2],1)),xyz],axis=-1)
    bvec=p-T
    D1=s0+2*dot(p,K);D2=s0-2*dot(bvec,K)
    if min(abs(D1).min(),abs(D2).min())<gap*(1-2e-12):
        raise ArithmeticError('Numerical node violated analytic off-root separation')
    Y,Z=slash(p+K)+M*I,slash(bvec-K)+m*I
    def V_(e):
        S=slash(e)
        return Y@S/D1[...,None,None]+S@Z/D2[...,None,None]
    def amp(e):return norm_spinors(p,bvec,M,m,V_(e))
    rows=(np.eye(4)[1:3] if name=='T' else np.eye(4)[:1])
    pnorm=sum(amp(e) for e in rows)
    ward=np.max(amp(K/k))/np.max(pnorm)
    if ward>2e-14 or np.min(pnorm)<=0.:
        raise ArithmeticError(('Coherent off-root current',ward))
    gauge=None
    if name=='L':
        alt=amp(np.array([0.,0.,0.,-k/w]))
        gauge=np.sum(abs(alt-pnorm))/np.sum(pnorm)
        if gauge>3e-10:raise ArithmeticError(('Spacelike equivalent density coupling',gauge))
    weighted=W*pnorm.mean(axis=1)/V
    A=weighted@np.repeat(f*(1-g),n)
    B=weighted@np.repeat(g*(1-f),n)
    D=weighted@np.repeat(f-g,n)
    if abs(A-B-D)>3e-14*max(A,B) or min(A,B)<=0.:
        raise ArithmeticError('Positive directional angular normalization')
    T0=.8
    pf,bf=(expit((.12-U)/T0) for U in (E,E-t))
    At,Bt=weighted@(pf*(1-bf)),weighted@(bf*(1-pf))
    therm=abs(Bt/At/np.exp(t/T0)-1)
    KMS_Bose=abs((At/Bt)*np.exp((q[0]-w)/T0)-1)
    ng=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],np.linalg.norm(q[1:]))-1]
    low,high,rho,noise=spectral(ch,k,w,name)
    net=(D+A/ng)*low-B*rho
    return dict(k=k,omega=w,cosine=c,mode=name,parent_energy_interval=limits,min_proven_denominator=gap,
                rms_ward_squared=ward,density_equivalence_error=gauge,thermal_pauli_directional_error=therm,full_bose_kms_balance_error=KMS_Bose,
                raw_gain_inverse_signed=[A,B,D],rho=rho,effective_noise_energy=noise,
                lesser_over_greater=low/high,raw_off_root_estimate_divided_by_n_scalar=net)


def run():
    out=[]
    for g in ('fine','third'):
        d,ch,_=load(g)
        for k,w,c in ((.0008,.0002,0.),(.00128,.00064,-.5),(.00256,.00064,.5)):
            for name in ('T','L'):
                for bounds in ((1.18,1.26),(1.85,2.25)):
                    a,b=(panel(d,ch,k,w,c,name,bounds,*z) for z in ((8,24),(16,56)))
                    err=float(max(abs(np.array(a['raw_gain_inverse_signed'])/
                                      np.array(b['raw_gain_inverse_signed'])-1.)))
                    if err>5e-9 or max(b['thermal_pauli_directional_error'],b['full_bose_kms_balance_error'])>2e-12:
                        raise ArithmeticError(('Off-root phase refinement',g,k,name,bounds,err,b))
                    b.update(dict(grid=g,mesh_error=err))
                    out.append(b)
    return dict(scope='Bounded spacelike leading collective model and bare off-root kernels; not matched or integrated.',data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
