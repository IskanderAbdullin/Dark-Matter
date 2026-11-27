"""Conserving versus loss-only toy line response under angular redraws."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre


def response(d, nu, eps=1e-10, spread=1.):
    z = eps + nu - 1j*d
    if spread == 0:
        return 1/z, 1/(eps-1j*d)
    a = (np.log(z + 1j*spread)-np.log(z - 1j*spread))/(2j*spread)
    return a, a/(1-nu*a)


def legendre_line(d, nu, lmax, smooth, eps=1e-10):
    # Unit half-span. Multiplication by mu couples neighboring orthonormal
    # Legendre functions. Schur complements form a stable continued fraction.
    a = np.arange(1., lmax+1.)
    c = a*a/(4*a*a-1)
    damping = nu*a*(a+1)/2 if smooth else nu*np.ones_like(a)
    base = eps-1j*d
    v = base+damping[-1]
    for k in range(lmax-1, 0, -1):
        v = base+damping[k-1]+c[k]/v
    return 1/(base+c[0]/v)


def check():
    z, w = roots_legendre(160)
    uniform = np.ones(z.size)
    dchecks = []
    for nu in (.2, 1., 3., 10.):
        for d in (0., .23, 1.5):
            diagonal = 0.001 + nu - 1j*(d-z)
            # Solve the full finite-dimensional gain/loss linear operator
            # independently rather than apply the analytic resolvent formula.
            matrix = np.diag(diagonal)-nu*np.outer(uniform, .5*w)
            solved = np.linalg.solve(matrix, uniform)
            expected = response(d,nu,.001)[1]
            dchecks.append(abs(.5*(w@solved)/expected-1.))
    assert max(dchecks) < 8e-12
    # A uniform probe is unaffected by a flavor-blind angular reorientation
    # when the detuning is angular independent.
    k = [response(.09,n,.005,0.)[1] for n in (.1,1.,30.)]
    assert max(abs(x/k[0]-1) for x in k) < 1e-14
    return max(dchecks)


def run():
    defect = check()
    out = []
    angular_refinement = []
    for n in (.1,.3,2/np.pi,1.,3.,10.,30.,100.):
        vals = []
        for i in (0,1):
            def f(d):
                return response(d,n)[i].real
            half = brentq(lambda d:f(d)-.5*f(0.),1e-9,1e3)
            # Split around the collisionless edge for the dilute controls.
            norm = 2*(quad(f,0.,1.,epsabs=1e-8)[0]+quad(f,1.,np.inf,epsabs=1e-8)[0])/np.pi
            assert abs(norm-1.)<1e-7
            vals.append(dict(peak=f(0.),half_width=half,area_over_pi=norm))
        pole = n-1/np.tan(1/n) if n>2/np.pi else None
        if n>=30.:
            assert abs(vals[1]['half_width']*3*n-1.)<1e-3
        smooth = lambda d,m: legendre_line(d,n,m,True).real
        half = brentq(lambda d:smooth(d,48)-.5*smooth(0.,48),1e-9,1e3)
        hi = brentq(lambda d:smooth(d,96)-.5*smooth(0.,96),1e-9,1e3)
        assert abs(half/hi-1)<1e-9
        md = max(abs(legendre_line(d,n,160,False)/response(d,n)[1]-1.)
                 for d in (0.,.2,1.,1.5))
        assert md<1e-9
        angular_refinement.append(md)
        out.append(dict(nu_over_angular_halfspan=n, loss_only=vals[0],
                        rotational_diffusion_half_width=hi,
                        conserving_redraw=vals[1], isolated_decay_pole=pole))
    return dict(scope="fixed radial momentum; linear angular detuning; "
                "uniform Markov redraw with gain equals loss on the isotropic component; "
                "no physical gauge collision calibration, no radial diffusion, no KMS rate",
                finite_matrix_max_relative_error=defect,
                independent_legendre_max_relative_error=max(angular_refinement), rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
