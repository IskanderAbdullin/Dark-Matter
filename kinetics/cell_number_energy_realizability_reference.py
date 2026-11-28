"""Number/kinetic-energy realizability for one massive fermion shell."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog


def kinetic(q,m):
    return q*q/(np.hypot(q,m)+m) if q>0 else 0.


def moment(a,b,m):
    return (b-a)*quad(lambda t: (a+(b-a)*t)**2*kinetic(a+(b-a)*t,m),
                     0.,1.,epsabs=0.,epsrel=2e-12)[0]


def fd_bounds(a,b,m,fraction):
    """Kinetic moment limits at fixed occupancy fraction, common factors omitted."""
    if not (0 <= a < b and m >= 0 and 0 <= fraction <= 1):
        raise ValueError("unphysical shell or occupancy fraction")
    if fraction == 0:
        return (0.,0.)
    if fraction == 1:
        e=moment(a,b,m)
        return e,e
    volume=(b-a)*(a*a+a*b+b*b)
    lower_cut=np.cbrt(a*a*a+fraction*volume)
    upper_cut=np.cbrt(b*b*b-fraction*volume)
    return moment(a,lower_cut,m),moment(upper_cut,b,m)


def run():
    rng=np.random.default_rng(94726)
    cases=[]
    for (n,a,b,m) in [(n,*params) for n in (521,4096) for params in
                       [(0.,1e-6,.1),(1e-4,.13,1.),(.05,5.,.3),(0.,2.,0.)]]:
        qs=np.cbrt(a**3+(b**3-a**3)*((np.arange(n)+.5)/n))
        k=np.array([kinetic(q,m) for q in qs])
        perm=rng.permutation(n)
        scale=k.max()
        order=k[perm]/scale
        vol=(b-a)*(a*a+a*b+b*b)/3
        for f in (.03,.11,.47,.8,.97):
            lo,hi=fd_bounds(a,b,m,f)
            m0=linprog(order,A_eq=np.ones((1,n)),b_eq=[n*f],bounds=(0,1),method="highs")
            m1=linprog(-order,A_eq=np.ones((1,n)),b_eq=[n*f],bounds=(0,1),method="highs")
            assert m0.success and m1.success
            err=max(abs(m0.fun*scale*vol/n-lo),abs(-m1.fun*scale*vol/n-hi))/hi
            assert err<(1e-3 if n==521 else 5e-5)
            full=moment(a,b,m)
            flat=f*full
            assert lo<=flat<=hi
            cases.append(dict(microcells=n,qlo=a,qhi=b,mass=m,fill=f,kinetic_min=lo,kinetic_max=hi,
                              flat=flat,relative_discrete_LP_error=err))
    # A nonnegative Bose distribution at this same N can approach N*Kmax,
    # whereas finite fermion capacity forces occupancy over a finite region.
    lo,hi=fd_bounds(0.,1e-6,.1,.8)
    bose_hi=(1e-6)**3/3*.8*kinetic(1e-6,.1)
    assert hi<bose_hi
    maxerr={n:max(x["relative_discrete_LP_error"] for x in cases if x["microcells"]==n) for n in (521,4096)}
    assert maxerr[4096]<.2*maxerr[521]
    return dict(cases=cases,maximum_relative_discrete_LP_error_by_grid=maxerr,
                fermion_to_unbounded_kinetic_max_at_near_rest_fill_p8=hi/bose_hi)

if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
