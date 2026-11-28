"""Sequential positivity-limited pair projection stress test.

A common event limiter protects both daughters' available populations and
the fermion's number capacity. It audits conservation axioms, not rate accuracy.
"""
import json
import math
import numpy as np
from pair_projection_reference import bracket_energy


def step(nb,nf,cap,ib,jb,wb,if_,jf,wf,signed_request):
    """Apply one signed count request to daughter populations in place."""
    ibs, fbs, wb2, wf2 = (ib,jb), (if_,jf), (1-wb,wb), (1-wf,wf)
    a = abs(signed_request)
    positive = signed_request >= 0
    if not positive:
        for ix, w in zip(ibs,wb2):
            if w > 0:
                a = min(a, max(0.,nb[ix])/w)
    for ix, w in zip(fbs,wf2):
        if w > 0:
            available = cap[ix]-nf[ix] if positive else nf[ix]
            a = min(a, max(0.,available)/w)
    rate = a if positive else -a
    for ix,w in zip(ibs,wb2):
        nb[ix] += rate*w
    for ix,w in zip(fbs,wf2):
        nf[ix] += rate*w
    return rate


def run():
    seed, n = 273426, 20000
    rng = np.random.default_rng(seed)
    eb, ef = np.geomspace(.01,6,96),np.geomspace(.02,7,96)
    a,b=np.exp(rng.uniform(np.log(eb[0]),np.log(eb[-1]),n)),np.exp(rng.uniform(np.log(ef[0]),np.log(ef[-1]),n))
    modes = {
        "increasing":np.full(n,.2),
        "reverse":np.full(n,-.2),
        "mixed":.2*rng.choice([-1.,1.],n),
        "zero_net":np.zeros(n)}
    ijw = (*bracket_energy(a,eb),*bracket_energy(b,ef))
    out={}
    for label, requests in modes.items():
        cap=np.geomspace(.01,1,96)
        nb, nf = .4*cap,.65*cap
        startb,startf=nb.copy(),nf.copy()
        peak=0.
        accepted=np.empty(n)
        for t in range(n):
            accepted[t]=step(nb,nf,cap,*(v[t] for v in ijw),requests[t])
            peak=max(peak,-float(nb.min()),-float(nf.min()),float((nf-cap).max()))
        assert np.all(accepted*requests >= 0) and np.all(abs(accepted) <= abs(requests))
        if label=="mixed":
            revb, revf = startb.copy(), startf.copy()
            for t in range(n-1,-1,-1):
                step(revb,revf,cap,*(v[t] for v in ijw),requests[t])
            order_dependence = math.fsum(abs(nf-revf))/math.fsum(cap)
        dn=math.fsum(accepted)
        denergy=math.fsum(accepted*(a+b))
        dn_error=max(abs(math.fsum(nb-startb)-dn),abs(math.fsum(nf-startf)-dn))
        e_error=abs(math.fsum((nb-startb)*eb)+math.fsum((nf-startf)*ef)-denergy)
        norm=math.fsum(abs(accepted))
        assert peak < 2e-15 and dn_error <= 2e-13*max(1.,norm)
        assert e_error <= 2e-13*max(1.,math.fsum(abs(accepted*(a+b))))
        if label=="zero_net":
            assert e_error==0 and dn_error==0 and np.all(nb==startb) and np.all(nf==startf)
        else:
            assert np.any(abs(accepted)<abs(requests))
        out[label]=dict(limited_events=int(np.count_nonzero(abs(accepted)<abs(requests))),
                        max_absolute_bound_violation=peak,absolute_number_residual=dn_error,
                        absolute_energy_residual=e_error,
                        accepted_absolute_fraction=norm/max(1.,math.fsum(abs(requests))))
    return dict(seed=seed,events_per_case=n,statistics=out,
                mixed_order_fermion_L1_over_total_capacity=order_dependence)

if __name__ == "__main__":
    print(json.dumps(run(),indent=2))
