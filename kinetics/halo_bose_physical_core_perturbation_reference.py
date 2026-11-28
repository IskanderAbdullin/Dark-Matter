"""Stable small core-difference integration in one common physical radius."""
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar
from scipy.special import gamma, zeta, roots_legendre
from halo_bose_primordial_cap_core_reference import integrate

G=gamma(1.5)
Ic=G*zeta(1.5)


def functions():
    # Outside the inner IR region, direct quadrature supplies I and I'.
    x=np.geomspace(.18,24.,2400)
    f=PchipInterpolator(x,np.array([integrate(y,0.) for y in x])/Ic)
    p=PchipInterpolator(x,np.array([-.5*integrate(y,0.,-.5) for y in x])/Ic)
    # For y<.2 analytic Taylor coefficients plus
    # the singular square-root term suffice in the derivative equation.
    # scipy's zeta requires explicit negative values from mpmath.
    from mpmath import zeta as mz
    zz=[float(mz(1.5-k))*(-1)**k for k in range(10)]
    import math
    zz=np.array([a/math.factorial(k) for k,a in enumerate(zz)])
    prim=np.polynomial.polynomial.polyder(zz)
    for y in (.19,.2):
        assert abs(f(y)-(G/Ic)*(np.polynomial.polynomial.polyval(y,zz)-2*np.sqrt(np.pi*y)))<3e-9
        assert abs(p(y)-(G/Ic)*(np.polynomial.polynomial.polyval(y,prim)-np.sqrt(np.pi/y)))<3e-8
    def fun(y):
        if y<.2:
            return G/Ic*(np.polynomial.polynomial.polyval(y,zz)-2*np.sqrt(np.pi*y))
        return float(f(y))
    def der(y):
        if y<.2:
            return G/Ic*(np.polynomial.polynomial.polyval(y,prim)-np.sqrt(np.pi/y))
        return float(p(y))
    return fun,der


def damping(y,d):
    if y>=d:
        return 0.
    # Exact cap deficit, analytic subtraction of Bose terms avoids
    # differencing two order-unity densities as F becomes large.
    rat=y/d
    top=np.sqrt(1-rat)
    def func(u):
        w=rat+u*u
        # d*(B(dw)-B(d)) is well behaved after multiplication by u²
        return 2*u*u*d*np.exp(d*w)*np.expm1(d*(1-w))/(np.expm1(d*w)*np.expm1(d))
    return np.sqrt(d)*quad(func,0.,top,epsabs=2e-11,epsrel=2e-10)[0]/Ic


def correction(w):
    if w<=0:
        return 4/3
    if w>=1:
        return 0.
    a=np.sqrt(1-w)
    return 2*a-2*np.sqrt(w)*np.arctan(a/np.sqrt(w))-2*a**3/3


def inner_speed(w):
    J=quad(lambda u:u*u*correction(u*u),0.,min(1.,w),epsabs=1e-12)[0]
    return -3*J/(np.sqrt(2)*Ic*w*w)


def run():
    F,Fp=functions()
    tiny=1e-10
    test=np.exp(2*np.log(tiny))
    log=(-np.log(1/tiny),np.log(16.))
    def base(l,z):
        return z[1],np.exp(2*l)*F(z[0])-z[1]
    ref=solve_ivp(base,log,[test/6,test/3],rtol=5e-12,atol=5e-15,
                  dense_output=True,max_step=.018)
    assert ref.success
    # Inner-to-outer point forcing amplitude, common physical scale.
    A=8/(15*Ic)
    k=G*np.sqrt(6*np.pi)/Ic
    def tangent(l,u):
        y=ref.sol(l)[0]
        return [u[1],np.exp(2*l)*Fp(y)*u[0]-u[1]]
    outer=solve_ivp(tangent,log,[-A, A*k*tiny/2],rtol=1e-11,atol=1e-13,dense_output=True)
    grid=np.array([.12,.337,1.,2.,4.,8.,12.])
    cases=[]
    im=minimize_scalar(inner_speed,bounds=(.03,1.),method='bounded')
    assert im.success
    root,w=roots_legendre(5)
    for cap in (100.,7000.,200000.,1e7,1e9):
        d=np.log1p(1/cap)
        a=d**1.5
        v=-damping(0.,d)/a
        def delta(l,z):
            y=ref.sol(l)[0]
            dz=z[0]*a
            yp=y+dz
            assert min(y,yp)>0
            mean=sum(ww*Fp(y+dz*(tt+1)/2)/2 for tt,ww in zip(root,w))
            r2=np.exp(2*l)
            return z[1], r2*(mean*z[0]-damping(yp,d)/a)-z[1]
        rr=solve_ivp(delta,log,[test*v/6,test*v/3],rtol=8e-10,atol=1e-11,
                     max_step=.024,dense_output=True)
        assert rr.success
        dc=rr.sol(np.log(grid))[1]
        vc=np.sqrt(ref.sol(np.log(grid))[1])
        dv=a*dc/(np.sqrt(vc*vc+a*dc)+vc)
        pred=outer.sol(np.log(grid))[1]/(2*vc)
        def differ(w):
            x=np.sqrt(6*d)*w
            lc=np.log(x)
            base=ref.sol(lc)[1]
            dif=a*rr.sol(lc)[1]
            return dif/(np.sqrt(base)+np.sqrt(base+dif))
        inn=minimize_scalar(differ,bounds=(.03,1.),method='bounded',
                            options=dict(xatol=1e-10))
        assert inn.success
        cases.append(dict(cap=cap,common_radius_units=grid.tolist(),
                          exact_inner_peak_radius_relative_to_rounding_scale=inn.x,
                          exact_inner_peak_velocity_difference_over_V0=inn.fun,
                          exact_inner_peak_velocity_difference_for_F583_4_V0_km_s=44.1*inn.fun,
                          relative_to_V0_speed_shift=dv.tolist(),
                          scaled_delta_minus_one_point_five=(dv/a).tolist(),
                          relative_defect_from_asymptotic=float(np.linalg.norm(dv/a-pred)/np.linalg.norm(pred))))
    assert cases[-1]['relative_defect_from_asymptotic']<.012
    assert abs(cases[-1]['exact_inner_peak_velocity_difference_over_V0']/(np.log1p(1/1e9)*im.fun)-1)<.005
    return dict(asymptotic_inner_peak_radius_in_rounding_units=im.x,
                asymptotic_peak_inner_velocity_difference_over_V0_delta=im.fun,
                scope="Critical-versus-capped isothermal family at fixed physical mass and dispersion, "
                      "small-signal perturbation in a common physical scale; no fit or tidal truncation.",
                potential_offset_leading_coefficient=-A,common_radius_units=grid.tolist(),
                asymptotic_speed_response_scaled_by_delta_one_point_five=pred.tolist(),rows=cases)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
