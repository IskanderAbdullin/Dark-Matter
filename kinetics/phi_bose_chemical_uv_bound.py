"""Certified thermal forward Yukawa number-source UV tail bound.

The positive-integrand upper bound is independent of angular quadrature.
Only common-temperature T and positive constant Delta mu are covered.
"""
import math
from scipy.integrate import quad
from kinetics.phi_bose_chemical_weak_reference import M1, MPH, G1, width, analytic, kinematics
from scipy.special import expit
import numpy as np


def forward_tail_upper(pmin, temp, mu1, mu2, muphi):
    """Strict bound on integral of R_phi over parent p>=pmin, for Delta mu>0.

    Uses f1 <= exp((mu1-E1)/T), 1-f2<=1, 1+b<=1+b(mphi),
    m1/E1 * p**2 <= m1*p, and E1>=p. Natural GeV units.
    """
    assert pmin>=0 and temp>0 and muphi<MPH
    delta=mu1-mu2-muphi
    assert delta>0
    bose_max=1/math.expm1((MPH-muphi)/temp)
    c=G1*width()*M1/(2*math.pi**2)
    return (c * math.exp(mu1/temp) *
            (1+bose_max)*(-math.expm1(-delta/temp)) *
            temp*(pmin+temp)*math.exp(-pmin/temp))


def finite_tail_scan(pmaxes=(2.,8.,12.,16.,20.,24.),temp=1.,
                     mu1=-.09,mu2=.05,muphi=-.2):
    """Finite-tail direct subtraction with certified *infinite* remainder."""
    assert temp==1., "current analytic() reference uses T1=T0=1 in this scan"
    pcut=60.
    full=analytic(pcut,temp,temp,mu1,mu2,muphi)['Rphi']
    rows=[]
    for pmax in pmaxes:
        assert 0<=pmax<pcut
        partial=analytic(pmax,temp,temp,mu1,mu2,muphi)['Rphi']
        # Difference of two independently integrated positive intervals is
        # not used as a high-precision tail when the tail is tiny.
        from kinetics.phi_bose_chemical_weak_reference import kinematics, moments
        def integrand(p):
            e,lo,hi,b=kinematics(p)
            f1=expit(-(e-mu1)/temp)
            aff=(mu1-mu2-muphi)/temp
            if b<1e-14:
                inv=expit(-((e-lo)-mu2)/temp)/math.expm1((lo-muphi)/temp)
                i0=inv
            else:
                i0=moments(e,lo,hi,temp,mu2,muphi)[0]/(2*b)
            return (G1*p*p/(2*math.pi**2)*width()*M1/e*
                    (1-f1)*math.expm1(aff)*i0)
        tail_to_60=quad(integrand,pmax,pcut,epsabs=1e-28,
                        epsrel=1e-9,limit=300)[0]
        bound=forward_tail_upper(pmax,temp,mu1,mu2,muphi)
        remainder=forward_tail_upper(pcut,temp,mu1,mu2,muphi)
        rows.append(dict(pmax=pmax,rate=partial,tail_to_60=tail_to_60,
                         strict_tail_bound=bound,
                         upper_tail_fraction=bound/partial,
                         finite_60_remainder_bound=remainder,
                         relative_finite_60_tail=tail_to_60/full))
    return rows


if __name__=="__main__":
    import json
    print(json.dumps(finite_tail_scan(),indent=2))
