"""Finite occupation regularization of a conditional isothermal critical Bose halo."""
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma, zeta

NCRIT=gamma(1.5)*zeta(1.5)

def occ(t,d):
    x=max(t,d)
    return np.exp(-x)/(-np.expm1(-x))

def integrate(y,d,power=.5):
    # For power=-1/2 remove the endpoint by t=u^2.
    if power==-.5:
        return 2*quad(lambda u:occ(u*u+y,d),0.,9.,points=[np.sqrt(max(d-y,0.))],
                      epsabs=2e-10,epsrel=5e-11)[0]
    return quad(lambda u:2*u*u*occ(u*u+y,d),0.,9.,
                points=[np.sqrt(max(d-y,0.))],epsabs=2e-11,epsrel=5e-11)[0]

def universal(w):
    if w>=1:
        return np.pi*np.sqrt(w)-4/3
    q=np.sqrt(1-w)
    return 2*np.sqrt(w)*np.arcsin(np.sqrt(w))+2*q-2*q*q*q/3-4/3

def profile(d,n0,high=False):
    yg=np.r_[0.,np.geomspace(d*1e-6,6.,1250 if high else 720)]
    fn=PchipInterpolator(yg,[integrate(y,d)/n0 for y in yg])
    eps=min(1e-6,np.sqrt(d)*1e-3)
    def hit(x,s):
        return s[0]-d
    out=solve_ivp(lambda x,s:[s[1]/x**2,x*x*fn(s[0])],(eps,2.2),
                  [eps**2/6,eps**3/3],events=hit,rtol=2e-9,atol=1e-13,dense_output=True)
    assert out.success
    vals=out.sol(np.array([.01,.03,.1,.3,1.,2.]))
    return out.t_events[0],fn(vals[0]),vals[1]/np.array([.01,.03,.1,.3,1.,2.])

def run():
    result=[]
    assert abs(integrate(0.,0.)/NCRIT-1)<2e-13
    for F in (1.,100.,5000.,7000.,200000.):
        d=np.log1p(1/F)
        n=integrate(0.,d)
        grad=-.5*integrate(0.,d,-.5)
        y=d*1e-5
        assert abs((integrate(y,d)-n)/(y*grad)-1)<1e-4
        a,b,c=profile(d,n);A,B,C=profile(d,n,True)
        assert np.max(abs(c/C-1))<4e-6 and np.max(abs(b/B-1))<4e-6
        if len(a):
            assert len(A)==1 and abs(a[0]/A[0]-1)<4e-6
        vals=[]
        for w in (.01,.1,1.,3.,10.):
            exact=(n-integrate(d*w,d))/np.sqrt(d)
            leading=universal(w)
            vals.append(dict(w=w,exact_scaled_density_drop=exact,
                             inner_scaling_limit=leading,relative=exact/leading))
        result.append(dict(cap=F,minimum_Bose_exponent=d,central_density_over_critical_at_fixed_V0=n/NCRIT,
                           density_quadratic_coefficient_at_fixed_central_scale=-grad/(6*n),
                           leading_rounding_radius_over_R=np.sqrt(6*d),
                           numerical_y_equals_delta_radius_over_R=float(A[0]) if len(A) else None,
                           matched_profile_radius_over_R=[.01,.03,.1,.3,1.,2.],
                           density_over_center=B.tolist(),velocity_over_V0_squared=C.tolist(),
                           matching_controls=vals))
    return dict(scope="Illustrative capped monotone energy-dependent scalar DF, spherical "
                      "self-gravitating isothermal local core, no finite halo fit, baryons, "
                      "mean-field interaction, collisionless formation or observational inference.",
                reference_critical_radial_linear_coefficient=np.pi/(np.sqrt(6)*NCRIT),
                cases=result)

if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
