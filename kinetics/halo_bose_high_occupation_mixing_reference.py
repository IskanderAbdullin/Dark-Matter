"""High-occupation excess-number requirements for local Bose halo cores."""
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit, gamma, roots_legendre
from halo_bose_primordial_cap_core_reference import integrate


def profile(F,kind):
    d=np.log1p(1/F)
    yg=np.r_[0.,np.geomspace(1e-12,2.,800),d]
    yg=np.unique(yg)
    ys=yg+d if kind=='shift' else yg
    den=PchipInterpolator(yg,[integrate(y,d) for y in ys])
    n=den(0.)
    eps=1e-7
    out=solve_ivp(lambda x,v:[v[1]/(x*x),x*x*den(v[0])/n],
                  (eps,1.),[eps*eps/6,eps**3/3],rtol=2e-10,atol=1e-14,
                  dense_output=True)
    assert out.success
    norm=n*out.y[1,-1]
    return d,out.sol,norm


def halo_excess(F,kind,eta,prof):
    d,solution,N=prof
    extra_x=lambda t: np.log1p(np.expm1(d)*t/(1.-t)/(1.+np.expm1(d)))
    Emax=extra_x(eta)+ (d if kind=='clip' else 0.)
    upper=(brentq(lambda x:solution(x)[0]-Emax,1e-7,1.)
           if solution(1.)[0]>Emax else 1.)
    r,w=roots_legendre(88)
    R=(r+1)*upper/2
    y=solution(R)[0]
    T=np.maximum(Emax-y,0.)
    plevel=np.minimum(np.maximum(d-y,0.),T) if kind=='clip' else np.zeros_like(T)
    p=(r[None,:]+1)/2
    t=plevel[:,None]+(T-plevel)[:,None]*p*p
    X=(t+y[:,None]-(d if kind=='clip' else 0.)).clip(0.)
    drop=(np.exp(d)*np.expm1(X))/np.expm1(d+X)
    integrand=np.sqrt(t)*np.maximum(eta-drop,0.)*p*w
    D=(T-plevel)*integrand.sum(axis=1)+2./3*eta*plevel**1.5
    return float(upper/2*((w*R*R)@D)*F/N)


def primordial(eta):
    result={'waterbag':eta}
    for n in (.5,1.,2.):
        # Use log-domain root for the very low thresholds.
        u=(-np.log1p(-eta))**(1/n)
        I=gamma(3/n)/n
        g=quad(lambda s:s*s*(eta+np.expm1(-(u*s)**n)),0.,1.,
               epsabs=1e-17,epsrel=2e-10)[0]
        result[f'exp_n{n:g}']=u**3*g/I
    eps=.045
    f0=expit(1/eps)
    a=np.exp(-1/eps)
    norm=quad(lambda u:u*u*expit((1-u)/eps)/f0,0.,4.,epsabs=1e-13)[0]
    u=eps*np.log1p(eta*(1+a)/(a*(1-eta)))
    g=quad(lambda s:s*s*(eta-a*np.expm1(u*s/eps)/(1+a*np.exp(u*s/eps))),0.,1.,
           epsabs=2e-17,epsrel=1e-10)[0]
    result['logistic_compact_class_shape']=u**3*g/norm
    return result


def run():
    results=[]
    for F in (100.,7000.):
        profiles={k:profile(F,k) for k in ('shift','clip')}
        rec=[]
        for eta in (.9,.5,.1,.01,1e-4,1e-6,1e-8,1e-10):
            src=primordial(eta)
            a,b=(halo_excess(F,k,eta,profiles[k]) for k in ('shift','clip'))
            assert a>0 and b>0 and min(src.values())>0
            rec.append(dict(cap=F,fractional_threshold_deficit=eta,
                            primordial_excess_over_input_count=src,
                            shifted_core_excess_over_core_count=a,
                            clipped_core_excess_over_core_count=b,
                            required_associated_mass_over_core_shifted={k:a/v for k,v in src.items()},
                            required_associated_mass_over_core_clipped={k:b/v for k,v in src.items()}))
        A=F*3*np.sqrt(6)*np.pi/(96*(F+1)**3*profiles['shift'][2])
        assert abs(rec[-1]['shifted_core_excess_over_core_count']/(A*1e-40)-1)<8e-5
        assert abs(rec[-1]['clipped_core_excess_over_core_count']/rec[-2]['clipped_core_excess_over_core_count']*100-1)<2e-4
        eta=1e-12
        src=primordial(eta)
        b=halo_excess(F,'clip',eta,profiles['clip'])
        a=A*eta**4
        rec.append(dict(cap=F,fractional_threshold_deficit=eta,
                        primordial_excess_over_input_count=src,
                        shifted_core_leading_asymptote_over_core_count=a,
                        clipped_core_excess_over_core_count=b,
                        required_associated_mass_over_core_clipped={k:b/v for k,v in src.items()}))
        results+=rec
    return dict(scope="Necessary collisionless phase-space mixing test on the core contribution only. "
                      "Initial distributions uniform in a finite Lagrangian patch, same scalar species, "
                      "no particle production or collisional rethermalization. No halo formation inference.",
                comparisons=results)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
