"""Vacuum two-diagram neutral-scalar photoproduction at finite unequal masses."""
import json

import numpy as np
from scipy.special import roots_legendre

I=np.eye(4,dtype=complex)
SIG=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
Z=np.zeros((2,2))
G=np.array([np.diag([1,1,-1,-1])]+[np.block([[Z,s],[-s,Z]]) for s in SIG])
MET=np.array([1,-1,-1,-1])


def dot(p,q):
    return np.sum(MET*np.asarray(p)*np.asarray(q),axis=-1)


def slash(p):
    return np.einsum('...j,jab->...ab',np.asarray(p)*MET,G)


def adj(u):
    return G[0]@np.swapaxes(u.conj(),-1,-2)@G[0]


def spinors(p,mass):
    W=np.sqrt(p[...,0]+mass)[...,None,None]
    U=W*np.broadcast_to(np.eye(2),(*p.shape[:-1],2,2))
    L=np.einsum('...j,jab->...ab',p[...,1:],SIG)/W
    return np.concatenate([U,L],axis=-2)


def norm_spinors(p,b,M,m,V):
    left=np.swapaxes(spinors(b,m).conj(),-1,-2)@G[0]
    return np.sum(abs(left@V@spinors(p,M))**2,axis=(-1,-2))

def amplitudes(p,k,b,M,m,pol):
    """Spin summed, before spin/polarization averaging or coupling factors."""
    pk=dot(p,k)[...,None,None]
    bk=dot(b,k)[...,None,None]
    P,B,K=map(slash,(p,b,k))
    ep=slash(pol)
    # Use the on-shell form and exact simple denominators to minimize
    # cancellation near collinearity. These are both same-charge diagrams.
    A=(2*dot(p,pol)[...,None,None]*I+K@ep)/(2*pk)
    D=(2*dot(b,pol)[...,None,None]*I-ep@K)/(-2*bk)
    def norm(V):
        z=np.trace((B+m*I)@V@(P+M*I)@adj(V),axis1=-2,axis2=-1)
        return z.real
    X=A+D
    return norm(X),norm(A),norm(D)


def invariant_total(s,cosa,M,m,mu):
    """Polarization and spin sum reduced to scalar on-shell invariants."""
    root=np.sqrt(s)
    a=(s-M*M)/2
    eb=(s+m*m-mu*mu)/(2*root)
    pb=np.sqrt(np.maximum(0.,(eb-m)*(eb+m)))
    c=a/root*(m*m/(eb+pb)+pb*(1-cosa))
    d=(M*M+m*m-mu*mu)/2+a-c
    return invariant_form(a,c,d,M,m)


def invariant_form(a,c,d,M,m):
    v=(1/a-1/c)/2
    H2=M*M/a**2+m*m/c**2-2*d/(a*c)
    X=c*M*M/a+a*m*m/c-2*d
    return -4*(d+M*m)*H2-8*v*X+16*v*v*a*c


def com(s,c,M,m,mu):
    t=np.sqrt(s)
    p=(s-M*M)/(2*t)
    r=np.sqrt((s-(m+mu)**2)*(s-(m-mu)**2))/(2*t)
    u=np.sqrt(np.maximum(0,1-c*c))
    fermion=np.stack([np.broadcast_to(np.hypot(p,M),c.shape),0*c,0*c,-p+0*c],axis=-1)
    photon=np.stack([p+0*c,0*c,0*c,p+0*c],axis=-1)
    out=np.stack([np.broadcast_to(np.hypot(r,m),c.shape),r*u,0*c,r*c],axis=-1)
    scalar=fermion+photon-out
    return fermion,photon,out,scalar,r/p


def trace_fixture(M,m,mu):
    rows=[]
    for s in (.0026,.004,.015,.1,1.,4.):
        if s<=max(M*M,(m+mu)**2):
            continue
        c=np.r_[-.9998,-.99,-.5,0.,.5,.99,.9998]
        p,k,b,q,kfactor=com(s,c,M,m,mu)
        assert np.max(abs(dot(q,q)-mu*mu))<2e-15
        a=np.array([amplitudes(p,k,b,M,m,u) for u in (np.array([0,1,0,0]),np.array([0,0,1,0]))])
        total=a[:,0].sum(axis=0)
        assert np.min(total)>0.
        scalar_error=np.max(abs(invariant_total(s,c,M,m,mu)-total)/total)
        assert scalar_error<2e-8,scalar_error
        cov=sum(v*amplitudes(p,k,b,M,m,u)[0] for v,u in zip(-MET,np.eye(4),strict=True))
        ward=amplitudes(p,k,b,M,m,k)[0]
        gauge=np.stack([np.broadcast_to(x,k.shape)+.31*k for x in ([0,1,0,0],[0,0,1,0])])
        ga=sum(amplitudes(p,k,b,M,m,u)[0] for u in gauge)
        cov_err=np.max(abs(cov-total)/total)
        gauge_err=np.max(abs(ga-total)/total)
        err=max(cov_err,gauge_err,np.max(abs(ward)/total))
        assert err<1e-9,err
        # Independent unsimplified numerator check on explicit spatial pols.
        forward=np.zeros_like(total)
        P,B,K=map(slash,(p,b,k))
        for u in (G[1],G[2]):
            Q=(P+K+M*I)@u/(s-M*M)+u@(B-K+m*I)/(-2*dot(b,k)[...,None,None])
            forward+=np.trace((B+m*I)@Q@(P+M*I)@adj(Q),axis1=-2,axis2=-1).real
        d=np.max(abs(forward-total)/total)
        full_ward=(P+K+M*I)@K/(s-M*M)+K@(B-K+m*I)/(-2*dot(b,k)[...,None,None])
        dward=norm_spinors(p,b,M,m,full_ward)
        assert max(d,np.max(abs(dward)/total))<1e-9, (d,dward/total)

        z,w=roots_legendre(40)
        pz,kz,bz,_,rat=com(s,.8*z,M,m,mu)
        physical=sum(amplitudes(pz,kz,bz,M,m,u)[0] for u in
                     (np.array([0,1,0,0]),np.array([0,0,1,0])))
        # dσ/dcos(θ) before g² y², with 1/4 for initial spin and γ.
        xs=.8*(w@physical)*rat/(128*np.pi*s)
        rows.append(dict(s=s,hard_central_cross_section_without_g2y2=xs,
                         min_and_max_total=[total.min(),total.max()],
                         complete_over_incoherent=(total/a[:,1:].sum(axis=(0,1))).tolist(),
                         ward_and_polarization_relative_error=err,
                         independent_scalar_invariant_relative_error=scalar_error,
                         independent_numerator_relative_error=d,
                         independent_ward_relative_error=np.max(abs(dward)/total)))
    return rows


def gate(M,m,mu):
    max_ET=np.hypot(5.,M)+4.
    q=1e-5
    # 2 K.q = s + m_phi² - m_out²; P_total cannot exceed E_total.
    ceiling=m*m-mu*mu+2*max_ET*(np.hypot(q,mu)+q)
    bound=(1+mu*mu-m*m)/(2*max_ET)
    lower_at_s1=(bound*bound-mu*mu)/(2*bound)
    max_lab_angle=[]
    for p,k in ((.45,.5),(1.4,.5),(1.4,1.),(1.4,4.),(5.,4.)):
        E=np.hypot(p,M)+k
        # Sharper bound using |P_total|<=p+k and s monotonic in cosθ.
        sm=m*m-mu*mu+2*(E*np.hypot(q,mu)+(p+k)*q)
        cosmin=(M*M + 2*k*np.hypot(p,M)-sm)/(2*p*k)
        ang=None if cosmin>=1 else float(np.arccos(np.clip(cosmin,-1,1)))
        max_lab_angle.append(dict(p=p,k=k,min_cos=float(cosmin),max_angle_rad=ang))
    assert ceiling<.02 and lower_at_s1>.027
    return dict(parent_mass=M,daughter_mass=m,scalar_mass=mu,maximum_s_for_cold_scalar=ceiling,
                minimum_scalar_q_if_s_at_least_one=lower_at_s1,
                local_incoming_angular_bounds=max_lab_angle)


def massless_fixture():
    angle=np.r_[-.9999,-.95,-.5,0.,.5,.95,.9999]
    errors=[]
    for s in (.007,1.,33.):
        p,k,b,_,_=com(s,angle,0.,0.,0.)
        v=sum(amplitudes(p,k,b,0.,0.,e)[0] for e in np.eye(4)[1:3])
        u=-s*(1-angle)/2
        t=-s-u
        expect=-4*t*t/(s*u)
        errors.append(float(np.max(abs(v-expect)/expect)))
    assert max(errors)<2e-7
    return max(errors)


def run():
    return dict(massless_spin_sum_relative_defect=massless_fixture(),scope='Bare real photon, scalar Yukawa and equal Abelian charges; two tree amplitudes. '
                      'No in-medium photon, virtual radiation, soft-hard matching or kinetic rate.',
                cases=[dict(**gate(M,m,mu),matrix_checks=trace_fixture(M,m,mu))
                       for M,m in ((.05,.003),(.003,.05)) for mu in (.0008720559389,.0008672414472)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
