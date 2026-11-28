"""Fully Bose-weighted held one-way encounter kernel on selected halo orbits."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from halo_bose_orbital_elastic_clock_reference import (
    background, trajectory, R0, V0, SIGMA_OVER_M, GYR, constants,
)


def average(a,b,xmin):
    # 1+f(x)=1/(1-e^-x). Product of the two factors has an elementary
    # polar integral, using 1-exp(-2a) constant over outgoing directions.
    lim=1/(-np.expm1(-a))**2
    sel=b>1e-10*a
    den=-np.expm1(-xmin)
    rat=np.exp(-xmin)*(-np.expm1(-2*b))/den
    k=np.log1p(rat)/np.where(sel,b,1.)
    return np.where(sel,(1+k)/(-np.expm1(-2*a)),lim)


def check_angles():
    rows=[]
    for a,b in ((4.,3.),(.2,.1),(.99,.98999),(.004,.0039999),(.005,1e-16)):
        q=quad(lambda c:.5/((-np.expm1(-a-b*c))*(-np.expm1(-a+b*c))),
               -1.,1.,epsabs=1e-10,epsrel=8e-11,limit=500)[0]
        p=average(a,b,a-b)
        rows.append(abs(q/p-1.))
    assert max(rows)<2e-9
    return max(rows)


def quad_u(n):
    x,w=roots_legendre(n)
    cut=np.r_[0.,np.geomspace(1e-7,13.,41)]
    lo,hi=cut[:-1,None],cut[1:,None]
    return (((lo+hi)/2+(hi-lo)/2*x).ravel(),((hi-lo)/2*w).ravel())


def quad_mu(n):
    x,w=roots_legendre(n)
    cut=np.r_[0.,np.geomspace(1e-4,1.,12)]
    lo,hi=cut[:-1,None],cut[1:,None]
    return (((lo+hi)/2+(hi-lo)/2*x).ravel(),((hi-lo)/2*w).ravel())


def local(d,y,v,n0,order):
    u,w=quad_u(order)
    x,wx=quad_mu(order-1)
    U=u*u; V=v*v
    P=U+V
    Q=(u[:,None]*v*x)**2
    pos=np.sqrt((u[:,None]-v)**2+2*u[:,None]*v*(1-x))
    neg=np.sqrt(P[:,None]+2*u[:,None]*v*x)
    S=np.sqrt(np.maximum(P[:,None]**2-4*Q,0.))
    a0=y+d
    a=a0+.25*P[:,None]
    zmin=a0+Q/(P[:,None]+S)
    bose=average(a,.25*S,zmin)
    f=np.exp(-a0-U/2)/(-np.expm1(-a0-U/2))
    p=w*U*f/(np.sqrt(2)*n0)
    num=p@(bose*(pos+neg)/2@wx)
    kick=p@(bose*(pos**3+neg**3)/4@wx)
    energy=(p*(U-V)/4)@(bose*(pos+neg)/2@wx)
    return np.array([num,kick,energy])


def orbit(pot,d,n0,rp,ra,N,M):
    if rp==ra:
        y,v2=pot(ra)
        return local(d,y,np.sqrt(v2),n0,M)
    yp,ya=pot([rp,ra])[0]
    L2=2*(ya-yp)/(rp**-2-ra**-2)
    E=ya+L2/(2*ra**2)
    x,w=roots_legendre(N)
    th=np.pi*(1+x)/2
    span=(ra-rp)/2
    r=(ra+rp)/2+span*np.cos(th)
    y=pot(r)[0]
    v2=2*(E-y)
    vr2=v2-L2/r**2
    q=w*np.sin(th)/np.sqrt(vr2)
    vals=np.array([local(d,a,b,n0,M) for a,b in zip(y,np.sqrt(v2),strict=True)])
    return q@vals/q.sum()


def run():
    c=check_angles()
    rho0=V0**2/(4*np.pi*constants.G*R0**2)
    factor=10*rho0*SIGMA_OVER_M*V0*GYR
    rows=[]
    for F in (100.,7000.):
        d=np.log1p(1/F)
        pot,n0=background(d)
        for rp,ra in ((.1,.1),(.3,.3),(1.,1.),(3.,3.),(10.,10.),(20.,20.),
                      (.1,3.),(.1,10.),(.1,20.),(5.,20.)):
            a=orbit(pot,d,n0,rp,ra,68,9)
            b=orbit(pot,d,n0,rp,ra,110,13)
            rel=np.max(abs(a-b)/(abs(b)+3e-5))
            assert rel<3e-3, (F,rp,ra,a,b,rel)
            vac=trajectory(pot,d,n0,rp,ra,100,14)['orbit_average_vacuum_moments']
            assert b[0]>vac[1] and b[1]>vac[2]
            rows.append(dict(peak_occupation=F,pericenter_R=rp,apocenter_R=ra,
                             vacuum_orbit_encounters_ten_Gyr=factor*vac[1],
                             stimulated_orbit_encounters_ten_Gyr=factor*b[0],
                             stimulated_to_vacuum_encounters=b[0]/vac[1],
                             stimulated_kick2_over_V0_squared_ten_Gyr=factor*b[1],
                             stimulated_tagged_energy_change_over_mV0_squared_ten_Gyr=factor*b[2],
                             joint_quadrature_max_relative_defect=rel))
    return dict(scope="Held thermal Bose bath, spherical potential, isotropic contact scattering; "
                      "one-way tagged collision and kick moments include both final Bose factors. "
                      "Not net gain, relaxation eigenvalue, or time evolution of tagged orbits.",
                exact_outgoing_polar_control_max_relative_defect=c,results=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
