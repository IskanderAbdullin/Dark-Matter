"""Sharp shell moment feasibility for a prescribed per-spin fermion cap."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from scalar_pair_expanding_continuation_reference import moments
from three_species_decay_onsager_reference import m
from warm_light_occupation_stage_audit_reference import digest, read

X, W = roots_legendre(32)
X, W = (X+1)/2, W/2


def volume(a, b):
    return (b-a)*(a*a+a*b+b*b)/3


def kinetic(p, mass):
    return p*p/(np.hypot(p, mass)+mass)


def extrema(a, b, mass, n, cap):
    """Kinetic limits for shell number n under 0 <= f <= cap."""
    if not (0 <= a < b and mass >= 0 and n >= 0 and 0 < cap <= 1):
        raise ValueError("invalid shell or cap")
    v = volume(a, b)
    if n > cap*v*(1+2e-14):
        return None
    if n == 0:
        return 0., 0.
    # 3 times the occupied volume. Avoid subtracting nearly equal roots:
    # (r-a)(r*r+r*a+a*a) == r**3-a**3.
    d = min(3*n/cap, 3*v)
    l = np.cbrt(a**3+d)
    h = np.cbrt(max(0., b**3-d))
    dl = d/(a*a+a*l+l*l)
    dh = d/(b*b+b*h+h*h)
    p, q = a+dl*X, b-dh*X
    return float(cap*dl*(W@(p*p*kinetic(p,mass)))), float(cap*dh*(W@(q*q*kinetic(q,mass))))


def feasible(a, b, mass, n, e, cap):
    lim = extrema(a, b, mass, n, cap)
    if lim is None:
        return False
    tol = 2e-11*max(abs(e), abs(lim[1]))
    return lim[0]-tol <= e <= lim[1]+tol


def minimum(a, b, mass, n, e):
    assert feasible(a, b, mass, n, e, 1.)
    if n == 0:
        return 0., 0.
    lo, hi = max(0., n/volume(a,b)*(1-4e-14)), 1.
    for _ in range(45):
        c = (lo+hi)/2
        if feasible(a,b,mass,n,e,c):
            hi = c
        else:
            lo = c
    return lo, hi


def independent_moments(a, b, mass, theta):
    p = a+(b-a)*X
    u = (kinetic(p,mass)-kinetic(a,mass))/(kinetic(b,mass)-kinetic(a,mass))
    w = (b-a)*W*p*p*expit(-theta[0]-theta[1]*u)
    return np.array([np.sum(w), w@kinetic(p,mass)])


def saturated_witness(a, b, mass, n, e):
    """Return inner cubed-volume fraction of a sliding unit-filled interval."""
    v = volume(a,b)
    # u is a fraction of the shell's total radial volume; the interval has
    # width n/v, so its number is n analytically regardless of this root.
    def mean(u):
        q = np.cbrt(a**3+3*v*u+3*n*X)
        return W@kinetic(q,mass)
    f = brentq(lambda u: mean(u)-e/n, 0., 1.-n/v, xtol=2e-15)
    return dict(inner_volume_fraction=f, interval_volume_fraction=n/v,
                f_essential_supremum=1., relative_kinetic_error=float(abs(mean(f)/(e/n)-1)))


def controls():
    n, e = .1, .09
    limits = extrema(0.,1.,0.,n,.5)
    assert np.isclose(limits[1], .5/4*(1-(1-3*n/.5)**(4/3)), rtol=3e-14)
    assert not feasible(0.,1.,0.,n,e,.5) and feasible(0.,1.,0.,n,e,1.)
    # The edge itself rounds to the original edge, but filled volume remains
    # finite and the kinetic moment divided by N approaches K(edge).
    n2=1e-27
    lo,hi=extrema(1.,1.0002,.02,n2,1.)
    assert abs(lo/n2/kinetic(1.,.02)-1)<3e-15
    assert abs(hi/n2/kinetic(1.0002,.02)-1)<3e-15
    assert not feasible(0.,1.,0.,.2,.15,.5)  # capacity failure
    return dict(massless_N=n, massless_E=e, massless_half_cap_Emax=limits[1],
                dilute_number=n2, dilute_kinetic_extrema=[lo,hi])


def run(directory):
    out=dict(controls=controls(), cases=[])
    for q in (16,32):
        past=read(directory/f'warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz')
        observed=read(directory/f'warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz')
        endpoints=past+[b["checkpoint"] for b in observed["blocks"]]
        results=[]
        for c in endpoints:
            s=c["snapshot"]
            assert c["snapshot_sha256"]==digest(s)
            p=np.asarray(s["lab_momentum_edges_GeV"])
            st=[s["entropy_variable_scalar"],s["entropy_variable_fermion"]]
            theta=np.asarray(st[1])
            ne=moments(p,st)[1].T*np.pi**2/2
            audit=np.asarray([independent_moments(a,b,m,t) for a,b,t in
                              zip(p[:-1],p[1:],theta.T,strict=True)])
            maxerr=float(np.max(np.abs(ne/np.maximum(audit,1e-280)-1)))
            assert maxerr < 2e-10
            caps=[]
            for a,b,(n,e) in zip(p[:-1],p[1:],ne,strict=True):
                assert feasible(a,b,m,n,e,.5)
                caps.append(minimum(a,b,m,n,e))
            ix=int(np.argmax(np.asarray(caps)[:,1]))
            peak=float(expit(-min(np.min(theta[0]),np.min(theta.sum(axis=0)))))
            assert caps[ix][1] < peak*(1+1e-7)
            witness=saturated_witness(p[ix],p[ix+1],m,*ne[ix])
            assert witness["relative_kinetic_error"]<2e-11
            results.append(dict(T_visible_GeV=s["T_visible_GeV"], snapshot_sha256=c["snapshot_sha256"],
                            maximum_of_shell_cap_brackets=np.max(caps,axis=0).tolist(), shell_index=ix,
                            reconstructed_peak=peak, independent_quadrature_max_relative_error=maxerr,
                            same_moment_saturated_shell=witness))
        out["cases"].append(dict(q=q,results=results))
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(parser.parse_args().tapes),indent=2,allow_nan=False))
