"""Static longitudinal HTL dipole kernel, not a moving-fermion linewidth."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import sici, roots_legendre

from early_pole_differential_phase_separation_reference import geometry


def kernel(x,avg):
    x = np.asarray(x)
    out = np.empty(x.shape)
    keep = abs(x)<.04
    y = x[keep]**2
    out[keep] = y*((1/18 if avg else 1/6)-y*((1/600 if avg else 1/120)
                       -y*(1/35280 if avg else 1/5040)))
    z = x[~keep]
    out[~keep] = 1-(sici(z)[0]/z if avg else np.sinc(z/np.pi))
    return out


def static_shape(y,avg,n,zmax=2e4):
    x,w=roots_legendre(n)
    # Smooth below the oscillation scale; above it refine to log intervals
    # small enough that oscillatory corrections are numerically immaterial.
    grid=np.r_[0.,np.geomspace(.01,zmax,220)]
    v=(grid[1:]+grid[:-1])[:,None]/2+np.diff(grid)[:,None]*x/2
    q=kernel(v*y,avg)*2*v/(1+v*v)**2
    estimate=np.sum(np.diff(grid)/2*(q@w))
    # For z>zmax the prefactor-integrated constant is exact. Bounding the
    # remaining oscillatory correction by |sin|<=1 or |Si|<=2 gives:
    tail=1/(1+zmax*zmax)
    error=(2 if avg else 1)*tail/(y*zmax)
    return estimate+tail,error


def run():
    rows=[]
    for g in (.001,.003,.004,.005,.01):
        dv,a=geometry(g)
        y=g*np.sqrt(2/3)*dv/a
        v,e=static_shape(y,True,80)
        nv,_=static_shape(y,True,160)
        check=quad(lambda u:static_shape(y*u,False,70)[0],0.,1.,epsabs=1e-7)[0]
        assert abs(v-nv)<2e-8 and abs(v-check)<3e-8
        rows.append(dict(g=g,separation_over_ideal_Debye_length=y,
                         static_phi=static_shape(y,False,160)[0],time_average_static_shape=nv,
                         integration_tail_error_bound=e,independent_time_average=check))
    # Only the static source is gauge-invariant via the rectangular loop.
    # Time-averaging at a linearly increasing length is a diagnostic functional,
    # not proof that a boosted moving dipole realizes that time dependence.
    return dict(scope='Laine et al. rectangular static potential shape hep-ph/0611300; '
                      'not the dynamical frequency kernel for ultrarelativistic lines.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
