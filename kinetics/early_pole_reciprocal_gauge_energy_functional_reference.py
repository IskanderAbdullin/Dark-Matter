"""Reciprocal first-order asymptotic gauge energy in the hard reference."""
import json
import numpy as np
from scipy.special import roots_legendre


def setup(n=96):
    u,w=roots_legendre(n)
    q=.25+(u+1)*7.75/2
    w=w*7.75/2
    f=np.array([1/np.expm1(q/.95+.4),1/(1+np.exp(q/.9+.3)),1/(1+np.exp(q/.45+1.))])
    # Degeneracy: two photon polarizations; four states per Dirac species.
    return q,w,f,np.r_[2.,4.,4.]


def state(a,m,g,q,w,f,d):
    p=q/a
    ds=w/a
    E=np.sqrt(p[None,:]**2+m[:,None]**2)
    J=np.sum((ds*p*p/E)*f,axis=1)
    N=(ds*p*p)[None,:]*d[:,None]/(2*np.pi**2)
    fac=g*g/np.pi**4
    # The photon and fermion derivatives are reciprocal by construction.
    I=fac*(.5*np.sum(J[1:]**2)+J[0]*np.sum(J[1:]))
    dJ=fac*np.r_[sum(J[1:]),J[0]+J[1:]]
    mass=2*np.pi**2*dJ/d
    disp=mass[:,None]/E
    trace=-dJ@(m*m*np.sum((ds*p*p/(E**3))*f,axis=1))
    return I,mass,disp,N,trace,J


def run(n):
    q,w,f,d=setup(n)
    g=.005
    out=[]
    for m in (np.zeros(3),np.array([0.,.05,.003])):
        I,A,H,N,R,J=state(1.,m,g,q,w,f,d)
        target=g*g/np.pi**2*np.r_[2*sum(J[1:]),J[0]+J[1:]]
        # A is one-half of the squared-mass increment.
        assert max(abs(2*A/target-1))<1e-13
        assert abs(np.sum(N*f*H)/(2*I)-1)<1e-13
        perturb=np.array([.2*np.sin(q),-.1*np.cos(2*q),.3*np.sin(1.5*q)])*f
        e=1e-5
        diff=(state(1.,m,g,q,w,f+e*perturb,d)[0]-state(1.,m,g,q,w,f-e*perturb,d)[0])/(2*e)
        deriv=np.sum(N*H*perturb)
        assert abs(diff/deriv-1)<1e-8
        ln_a=(state(np.exp(e),m,g,q,w,f,d)[0]-state(np.exp(-e),m,g,q,w,f,d)[0])/(2*e)
        trace_error=abs(ln_a+4*I-R)/I
        assert trace_error<2e-8
        out.append(dict(masses=m.tolist(),interaction_density=I,
                        photon_and_fermion_mass_sq=(2*A).tolist(),
                        trace_interaction=R,physical_scale_derivative_normalized=ln_a/I,
                        difference_trace_defect=trace_error,first_variation_ratio=diff/deriv))
    return dict(scope='First-order massless-loop energy functional, with optional reciprocal finite-bare-energy '
                      'shape proxy in a hard domain. Not a matched massive interaction free energy, collision '
                      'operator, or coupled expansion trajectory.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(96),indent=2,allow_nan=False))
