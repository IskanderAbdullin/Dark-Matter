"""Moment-preserving inversion at fixed reconstruction temperature (reference)."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, logit


def shell_average(q0, q1, mass, temp, shift, statistics, holes=False):
    """Average occupation, with shift=chemical potential minus rest mass."""
    d = q1-q0
    volume_density = (q1*q1+q0*q1+q0*q0)/3
    if statistics == "BE" and shift > 0:
        raise ValueError("positive rest-mass-relative BE chemical potential")
    def kinetic(q):
        return q*q/(np.hypot(q,mass)+mass) if q > 0 else 0.
    def integrand(t):
        q = q0+d*t
        val = (kinetic(q)-shift)/temp
        if statistics == "FD":
            return q*q/volume_density*expit(val if holes else -val)
        if q == 0 and val == 0:
            return 2*mass*temp/volume_density
        with np.errstate(over="ignore"):
            return q*q/volume_density/np.expm1(val)
    points = []
    if statistics == "BE" and shift < 0:
        transition = np.sqrt((-shift)*(-shift+2*mass))
        if q0 < transition < q1:
            points.append((transition-q0)/d)
    return quad(integrand,0.,1.,points=points,epsrel=1e-10,epsabs=0.)[0]


def infer_shift(q0, q1, mass, temp, target, statistics):
    if not (0 <= q0 < q1 and mass >= 0 and temp > 0 and np.isfinite(target)):
        raise ValueError("invalid shell parameters")
    if statistics not in ("FD","BE"):
        raise ValueError("unknown statistics")
    def kinetic(q):
        return q*q/(np.hypot(q,mass)+mass) if q > 0 else 0.
    a, b = kinetic(q0), kinetic(q1)
    if statistics == "FD":
        if not 0 < target < 1:
            raise ValueError("exact saturation/empty do not have finite FD shifts")
        lo, hi = a+temp*logit(target), b+temp*logit(target)
        holes = target > .5
        def f(s):
            avg = shell_average(q0,q1,mass,temp,s,statistics,holes)
            return (1-target)-avg if holes else avg-target
    else:
        if not target > 0:
            raise ValueError("exactly empty cells do not have a finite BE shift")
        hi = 0.
        limit = shell_average(q0,q1,mass,temp,hi,statistics)
        if target > limit:
            raise ValueError("occupancy exceeds this thermal BE cell's capacity")
        if target == limit:
            return hi
        lo = a-temp*np.log1p(1/target)
        def f(s):
            return shell_average(q0,q1,mass,temp,s,statistics)-target
    # Differences in sufficiently narrow cells may round below the mu scale.
    endpoint_tol = 1e-13*(min(target,1-target) if statistics == "FD" else target)
    for s in (lo,hi):
        if abs(f(s)) <= endpoint_tol:
            return s
    return brentq(f,lo,hi,xtol=np.nextafter(0.,1.),rtol=9e-16,maxiter=160)


def run():
    mass, temp = .0004, 1.
    rows = []
    for stat, shifts in [("BE",[-1.,-2e-6,-1e-12,0.]),
                         ("FD",[-26.,-1.,.1,14.])]:
        for shift in shifts:
            recovered, moments = [], []
            for lo,hi in [(0.,1e-7),(1e-7,4e-6),(4e-6,1e-4),(1e-4,.1),(.1,3.)]:
                n = shell_average(lo,hi,mass,temp,shift,stat)
                s = infer_shift(lo,hi,mass,temp,n,stat)
                recovered.append(abs(s-shift))
                moments.append(abs(shell_average(lo,hi,mass,temp,s,stat)/n-1))
            rows.append(dict(stat=stat,shift=shift,largest_abs_shift_error=max(recovered),
                             largest_relative_moment_residual=max(moments)))
    assert max(r["largest_relative_moment_residual"] for r in rows) < 1e-10
    assert max(r["largest_abs_shift_error"] for r in rows if r["stat"]=="BE") < 1e-10
    limit=shell_average(0.,1e-4,mass,temp,0.,"BE")
    for case in [(1.001*limit,"BE"),(1.,"FD"),(0.,"FD")]:
        try:
            infer_shift(0.,1e-4,mass,temp,*case)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid thermal moment accepted")
    from numpy.polynomial.legendre import leggauss
    nodes, weights = leggauss(160)
    low, high, shift_ref = 0., 1.e-7, -1.e-12
    q = low+(high-low)*(nodes+1)/2
    eps = q*q/(np.hypot(q,mass)+mass)
    moments = np.sum(weights*q*q/np.expm1((eps-shift_ref)/temp))/np.sum(weights*q*q)
    inferred = infer_shift(low,high,mass,temp,moments,"BE")
    independent_relative_error = abs(inferred/shift_ref-1)
    assert independent_relative_error < 1e-9
    return {"mass":mass, "independent_BE_near_rest_relative_shift_error":independent_relative_error,"temperature":temp,"rows":rows,"rejected_above_capacity":True}

if __name__ == "__main__":
    print(json.dumps(run(),indent=2))
