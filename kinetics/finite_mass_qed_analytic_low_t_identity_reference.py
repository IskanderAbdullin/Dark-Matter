"""Analytic leading isotropic low-T photon-only cancellation at finite masses."""
import json

import numpy as np
import sympy as sp

from finite_mass_qed_matched_low_temperature_photon_bath_reference import real, shifted_rate


def primitive(w,r,u,m,a):
    B=m*m
    return (-4*(u+a)*(w-B/w-2*u*np.log(w))+
            8*r*(u+a)*(w-np.log(w))+
            r*r*(2*w*w-8*w+4*np.log(w)))


def exact_real(m,s,mode,r):
    a=m if mode=='scalar' else 0.
    u=(1+m*m-s*s)/2
    p=np.sqrt((1-(m+s)**2)*(1-(m-s)**2))/2
    S=1-2*r
    h=(S+m*m-s*s+np.sqrt((S-(m+s)**2)*(S-(m-s)**2)))/(2*S)
    l=m*m/(S*h)
    return (primitive(h,r,u,m,a)-primitive(l,r,u,m,a))/(8*p*(u+a))


def proof():
    w,r,u,m,a,x=sp.symbols('w r u m a x', positive=True)
    b=m*m
    d=u+r*(w-1)
    H=1+b/w**2-2*d/w
    V=(-1+1/w)/2
    X=w+b/w-2*d
    U=sp.expand(-4*(d+a)*H-8*r*V*X+16*r*r*V*V*w)
    coefficients=[sp.expand(U).coeff(r,j) for j in range(3)]
    p_primitive=-4*(u+a)*(w-b/w-2*u*sp.log(w))+8*r*(u+a)*(w-sp.log(w))+r*r*(2*w*w-8*w+4*sp.log(w))
    assert sp.simplify(sp.diff(p_primitive,w)-U)==0
    # w solves S*w² - (S + m² - s²)*w + m² = 0. Differentiate at r=0.
    wp=w*(w-1)/(w-u)
    wpp=-(2*wp**2+2*(-4*w+2)*wp)/(2*(w-u))
    F=sp.integrate(2*coefficients[2],w)
    boundary=F+2*coefficients[1]*wp+sp.diff(coefficients[0],w)*wp**2+coefficients[0]*wpp
    y=m*m/x
    E=(x+y)/2
    p=(x-y)/2
    actual=(boundary.subs(w,x)-boundary.subs(w,y)).subs(u,E)/(8*p*(E+a))
    # h is the mass-shift parameter, dM_T/dh=2. The correction to
    # the parent spin numerator's time component has derivative 2 too.
    dE=4-2*E
    external=(E*dE-2)/p**2-4+(dE-sp.log(x/y)/p)/(E+a)
    # log(x/y) is the only independent logarithm. Both its coefficient
    # and the constant part of the sum reduce to zero.
    L=sp.Symbol('L',real=True)
    residual=sp.factor(sp.together(sp.expand(sp.expand_log(actual+external+2,force=True)
                                           .subs(sp.log(m),sp.log(x)-L/2))))
    if residual!=0:
        raise ArithmeticError(f'Analytic low-T identity failed: {residual}.')
    return str(residual),[str(c) for c in coefficients]


def run():
    residue,coefficients=proof()
    comparisons=[]
    for m,s in ((.003,0),(.03,.04),(.06,.0174),(.06,.4),(.2,.1),(.2,.65),(.5,.12)):
        for mode in ('left','scalar'):
            direct,B=real(m,s,mode,116)
            kmax=(1-(m+s)**2)/2
            points=[-.21,-.04,0.,.01,.45*kmax,.81*kmax]
            defects=[abs(exact_real(m,s,mode,r)-direct(r)) for r in points]
            if max(defects)>2e-9:
                raise ArithmeticError(f'Primitive failed the signed gamma-matrix real replay: {m,s,mode,defects}.')
            comparisons.append(dict(m_over_M=m,scalar_over_M=s,current=mode,
                                    maximum_signed_real_absolute_defect=float(max(defects)),
                                    points=points,external_K_and_phase_coefficient=shifted_rate(m,s,mode)['combined']))
    return dict(scope='Analytic O(T²) rest-frame photon-only equal-charge benchmark, for the scalar and chiral currents; away from mass/phase-space endpoints. Not a hot-medium rate.',
                symbolic_sum_plus_two=residue,integrand_r_power_coefficients=coefficients,
                direct_signed_real_checks=comparisons)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
