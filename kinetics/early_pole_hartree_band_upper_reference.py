"""Shape-independent conditional Hartree bound in a fixed primordial momentum band."""
import json
import numpy as np
from early_pole_decay_mode_clock_reference import G,X,KAPPA
from early_pole_shell_hartree_feedback_reference import S,ZGOAL

WMAX=27/16


def band(lam,p,fac,a,xend):
    # c is an assumed mass-to-band-edge ratio used to define a first
    # crossing of n0. Any larger c tightens the rate but raises n0.
    V=p**3/(6*np.pi**2)
    n=2*(a*p)**2*np.hypot(p,a*p)/lam
    R=np.hypot(1.,1/a)
    K=KAPPA*G*fac*R*(1+V/n)/lam
    top=n+K*(xend**5-X**5)/5
    when=(X**5+5*max(ZGOAL-n,0.)/K)**.2
    return top,n,when


def run():
    rr=np.linspace(0.,1.,2001)
    g=(1-rr)*(1+rr)**3
    assert abs(max(g)-WMAX)<1e-14
    rows=[]
    for l in (7.65175e-8,1.08515714e-6):
        for band_edge in (1e-5,1e-4):
            for b in (1.,WMAX):
                for xend in (.006,.008,.01,.02,.05):
                    # Any individual a is a valid upper, so take the tightest sampled one.
                    vals=[(band(l,band_edge,b,a,xend),a) for a in (.5,1.,2.,4.,8.,16.,32.)]
                    (n,n0,t),a=min(vals,key=lambda x:x[0][0])
                    rows.append(dict(lambda_self_quartic=l,band_q_upper=band_edge,
                                     parent_x_end=xend,Born_width_ratio_upper=b,
                                     first_crossing_mass_over_band_edge=a,
                                     first_crossing_number_over_T3=n0,
                                     number_band_upper_over_T3=n,yield_upper=n/S,
                                     fraction_upper_of_illustrative_target=n/ZGOAL,
                                     parent_x_necessary_to_reach_target_at_this_a=t))
    return dict(scope="Arbitrary scalar occupation shape below the fixed comoving momentum band, "
                      "adopted positive tadpole gap equation; radiation-clock FD parent, scalar-Yukawa "
                      "Born rest width, all scalar positive contributions retained and no other "
                      "source or redistribution into the band. Fermion daughter can only reduce "
                      "forward gain relative to being empty. No full abundance or field prediction.",
                uniform_width_upper_any_fermion_pole=WMAX,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
