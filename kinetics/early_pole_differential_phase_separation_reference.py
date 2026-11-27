"""Straight-path equal-charge Gaussian-noise phase-difference witness."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erf

import early_pole_finite_resonance_aperture_reference as z


def missing_corr(y, shape):
    # Dimensionless path-integrated 1-C, divided by total elapsed time.
    if shape == 'gaussian':
        if y<.02:
            return y*y/6-y**4/40+y**6/336-y**8/3456
        return 1-np.sqrt(np.pi/2)*erf(y/np.sqrt(2))/y
    if y<.001:
        return y/2-y*y/6+y**3/24-y**4/120
    return 1+np.expm1(-y)/y


def geometry(g):
    M, b = np.hypot(z.X,g/2), np.hypot(z.RB*z.X,g/2)
    E, l, h, _, _ = z.intervals(M,z.MS,b,z.Q)
    ep = (l+h)/2
    p = np.sqrt((ep-M)*(ep+M))
    eb = ep-E
    mu = (p*p+z.Q*z.Q+b*b-eb*eb)/(2*p*z.Q)
    assert abs(mu)<1
    v = np.array([0,0,p])/ep-np.array([-z.Q*np.sqrt(1-mu*mu),0,p-z.Q*mu])/eb
    a = z.Q*p/eb
    return np.linalg.norm(v), a


def run():
    out=[]
    err=0.
    for g in (.001,.003,.004,.005,.01):
        dv,a=geometry(g)
        # Only a length-scale comparison, using ideal massless two-Dirac Debye.
        eta = g*np.sqrt(2/3)*dv/a
        for shape in ('gaussian','exponential'):
            k=missing_corr(eta,shape)
            for y in np.geomspace(1e-5,20,12):
                f=(lambda u: u*u/2) if shape=='gaussian' else (lambda u:u)
                val=quad(lambda x:-np.expm1(-f(y*x)),0.,1.,epsabs=1e-14)[0]
                err=max(err,abs(missing_corr(y,shape)-val))
            for r in (.1,1.,10.):
                tc=brentq(lambda t:2*r*t*missing_corr(eta*t,shape)-1.,1e-8,1e5)
                out.append(dict(g=g,velocity_difference=dv,angular_halfspan=a,
                                separation_at_inverse_halfspan_over_length=eta, shape=shape,
                                imposed_single_phase_decay_exponent_per_inverse_halfspan=r,
                                pair_exponent_fraction_of_independent_at_t_halfspan=k,
                                pair_efold_time_over_independent=tc*2*r))
    assert err<2e-12
    return dict(scope='White-in-time scalar potential toy; not a dynamical gauge width',
                exact_equal_paths_pair_exponent=0.,independent_quad_max_absolute_error=err,rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
