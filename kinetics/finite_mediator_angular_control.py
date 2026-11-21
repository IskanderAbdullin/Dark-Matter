"""Propagator-only angular diagnostic, NOT a physical chi-SM rate.
Run: python kinetics/finite_mediator_angular_control.py
"""
import json
import math

def analytic(x):
    if x < 0:
        raise ValueError("x must be nonnegative")
    if x == 0:
        return 1., 1.
    isotropic = 1. / (1. + x)
    if x < 1e-3:
        transport = sum((-x)**n * 2*(n+1)/(n+2) for n in range(12))
    else:
        transport = 2*(math.log1p(x)-x/(1+x))/x**2
    return isotropic, transport

def simpson(x, n=32768):
    if n % 2:
        raise ValueError("even n required")
    h=1./n
    sums=[0.,0.]
    for i in range(n+1):
        u=i*h
        v=(1+x*u)**-2
        w=1 if i in (0,n) else (4 if i%2 else 2)
        sums[0]+=w*v
        sums[1]+=w*2*u*v
    return tuple(v*h/3 for v in sums)

def run():
    rows=[]
    for ratio in [0.1,0.3,1.,3.,8.59926345,10.,30.,100.]:
        x=4/ratio**2
        a=analytic(x)
        b=simpson(x)
        error=max(abs(i-j) for i,j in zip(a,b))
        assert error<2e-9, (ratio,error)
        assert 0<a[1]<=a[0]<=1
        rows.append({"mV_over_pcm":ratio,"x":x,
                     "isotropic_ratio":a[0],"transport_ratio":a[1],
                     "max_abs_quadrature_error":error,
                     "max_pointwise_contact_relative_defect":1-(1+x)**-2})
    assert abs(rows[4]["max_pointwise_contact_relative_defect"]-.1)<1e-9
    print(json.dumps({"scope":"propagator-only, constant CM angular numerator; not physical 5D thermal collision",
                      "quadrature":"composite Simpson 32768 subintervals",
                      "rows":rows},indent=2))

if __name__=="__main__":
    run()
