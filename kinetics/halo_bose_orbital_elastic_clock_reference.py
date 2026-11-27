"""Orbit-averaged elastic diagnostics in the fixed spherical Bose benchmark."""
import json
import numpy as np
from scipy import constants
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from halo_bose_primordial_cap_core_reference import integrate

R0=.6e3*constants.parsec
V0=44.1e3
SIGMA_OVER_M=.125  # 1.25 cm^2/g, illustrative source benchmark, in m^2/kg
GYR=1e9*365.25*24*3600


def background(d):
    ys=np.r_[0.,np.geomspace(1e-12,20.,1400),d]
    yg=np.unique(ys)
    # Use shifted Bose rather than a plateau; d=0 is the critical limit.
    den=np.array([integrate(y+d,d) for y in yg])
    f=PchipInterpolator(yg,den/den[0])
    l=-18.
    s=solve_ivp(lambda z,u:[u[1],np.exp(2*z)*f(u[0])-u[1]],(l,np.log(55.)),
                [np.exp(2*l)/6,np.exp(2*l)/3.],rtol=1e-10,atol=1e-13,
                max_step=.025,dense_output=True)
    assert s.success
    return lambda x:s.sol(np.log(np.maximum(x,np.exp(l)))), den[0]


def moments(d,y,v,n0,order):
    # t=u^2/2; retain log-segmented quadrature in target speed to resolve
    # the critical p=0 tail even within the first orbit's inner passage.
    t,w=roots_legendre(order)
    e=np.r_[0.,np.geomspace(1e-7,13.,36)]
    aa,bb=e[:-1,None],e[1:,None]
    u=((aa+bb)/2+(bb-aa)/2*t).ravel()
    du=((bb-aa)/2*w).ravel()
    def fold(y,v):
        z=.5*u*u+y+d
        f=np.exp(-z)/(-np.expm1(-z))
        q=du*u*u*f/(np.sqrt(2)*n0)
        big=np.maximum(u,v)
        small=np.minimum(u,v)
        mean=big+small*small/(3*big)
        mean3=big**3+2*big*small*small+small**4/(5*big)
        kin=(u*u-v*v)/4.
        return [q.sum(),q@mean,q@mean3/2,q@(mean*kin)]
    return np.asarray([fold(a,b) for a,b in zip(np.atleast_1d(y),np.atleast_1d(v),strict=True)])


def trajectory(pot,d,n0,rp,ra,order,morder):
    if rp==ra:
        y, v2=pot(ra)
        return dict(period_R_over_V0=None,mean_radius_R=ra,
                    orbit_average_vacuum_moments=moments(d,y,np.sqrt(v2),n0,morder)[0])
    if rp:
        yp, ya=pot([rp,ra])[0]
        L2=2*(ya-yp)/(rp**-2-ra**-2)
    else:
        L2=0.
    E=pot(ra)[0]+L2/(2*ra*ra)
    n,w=roots_legendre(order)
    th=np.pi*(1+n)/2
    mid=.5*(ra+rp);h=.5*(ra-rp)
    x=mid+h*np.cos(th)
    y=pot(x)[0]
    v2=2*(E-y)
    vr2=v2-L2/x**2
    assert np.min(vr2)>0.
    zz=w*np.sin(th)*h*np.pi/2/np.sqrt(vr2)
    obs=moments(d,y,np.sqrt(v2),n0,morder)
    return dict(period_R_over_V0=2*zz.sum(),mean_radius_R=float(zz@x/zz.sum()),
                orbit_average_vacuum_moments=(zz@obs)/zz.sum())


def run():
    rho0=V0*V0/(4*np.pi*constants.G*R0*R0)
    clock=rho0*SIGMA_OVER_M*V0*GYR
    halo=[]
    for d in (0.,np.log1p(1/7000)):
        pot,n0=background(d)
        cases=[]
        for a,b in ((.3,.3),(1.,1.),(3.,3.),(10.,10.),(20.,20.),
                    (.1,3.),(.1,10.),(1.,10.),(.1,20.),(1.,20.),(5.,20.),(0.,20.)):
            lo=trajectory(pot,d,n0,a,b,88,10)
            hi=trajectory(pot,d,n0,a,b,156,15)
            A=hi['orbit_average_vacuum_moments']
            B=lo['orbit_average_vacuum_moments']
            assert np.max(abs(A-B)/(1.e-6+abs(A)))<1e-3, (d,a,b,A,B)
            fp=dict(pericenter_R=a,apocenter_R=b,mean_radius_R=hi['mean_radius_R'],
                    period_Gyr=None if a==b else hi['period_R_over_V0']*R0/V0/GYR,
                    mean_density_over_central=A[0],expected_vacuum_tagged_encounters_in_ten_Gyr=10.*clock*A[1],
                    velocity_kick_squared_over_V0_squared_accumulated_in_ten_Gyr=10.*clock*A[2],
                    tagged_kinetic_energy_change_over_mass_V0_squared_ten_Gyr=10.*clock*A[3],
                    convergence_max_relative_rate_defect=float(np.max(abs(A-B)/(1e-6+abs(A)))))
            cases.append(fp)
        halo.append(dict(central_offset=d,orbits=cases))
    return dict(scope="Fixed spherical potential, illustrative F583-4 scale, constant isotropic vacuum elastic "
                      "cross section, uncorrelated fixed Bose bath. Encounters do not actually change the held "
                      "test orbit. Tagged classical moments are not kinetic relaxation eigenvalues or evidence "
                      "of thermalization. No Bose final-state factors, assembly or observational fit.",
                R_kpc=.6,V0_km_s=44.1,sigma_over_mass_cm2_per_g=1.25,
                central_density_g_cm3=rho0*1e-3,rate_frequency_per_Gyr=clock,profiles=halo)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
