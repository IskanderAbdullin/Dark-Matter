"""Quantify the inherited scalar occupation in the kinematically accessible band."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from scalar_pair_expanding_continuation_reference import moments
from three_species_decay_onsager_reference import m, mp
from warm_pair_selective_loss_bound_reference import sha


def band(p,b):
    p=float(p)/b
    E=np.hypot(m,p)
    c=m*m/(E+p)
    return b*np.array([np.sqrt((c-mp)*(c+mp)),np.sqrt((E+p-mp)*(E+p+mp))])


def number(q,th,lo,hi,n):
    a,b=q[:-1],q[1:]
    l=np.maximum(a,lo)
    h=np.minimum(b,hi)
    used=h>l
    j=np.flatnonzero(used)
    r,w=roots_legendre(n)
    x=((h[j]+l[j])[:,None]+(h[j]-l[j])[:,None]*r)/2.
    dx=(h[j]-l[j])[:,None]*w/2.
    K=lambda p:p*p/(np.hypot(p,mp)+mp)
    ka=K(a[j])
    kd=K(b[j])-ka
    z=th[0,j,None]+th[1,j,None]*(K(x)-ka[:,None])/kd[:,None]
    assert np.all(z>0.)
    scalar=np.exp(-z)/(-np.expm1(-z))
    ni=float(np.sum(dx*x*x*scalar)/(2*np.pi**2))
    if len(j):
        ends=th[0,j,None]+th[1,j,None]*(K(np.column_stack([l[j],h[j]]))-ka[:,None])/kd[:,None]
        mn=float(np.min(ends))
        limit=float(np.exp(-mn)/(-np.expm1(-mn)))
    else:
        limit=0.
    return ni,limit


def run(shapes):
    rows=[]
    for p in shapes:
        S=json.loads(p.read_text())["snapshots"][-1]
        q=np.asarray(S["lab_momentum_edges_GeV"])
        th=[np.asarray(S[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
        N=float(np.sum(moments(q,th)[0,0]))
        full,peak=number(q,th[0],q[0],q[-1],44)
        assert abs(full/N-1.)<2e-9
        data=[]
        lo_prev=0.
        for b in (1.,1.15,1.5,2.,3.,5.):
            lo,hi=band(q[-1],b)
            assert lo>=lo_prev and hi>=lo
            lo_prev=lo
            n,f=number(q,th[0],lo,hi,28)
            nn,ff=number(q,th[0],lo,hi,44)
            assert abs(n-nn)<max(3e-7*nn,3e-10*N) and abs(f-ff)<2e-14
            data.append(dict(b=b,outer_containing_lab_interval_GeV=(np.array([lo,hi])/b).tolist(),
                             outer_containing_anchor_q_interval_GeV=[lo,hi],
                             inherited_resolved_number_fraction_in_band=nn/N,
                             inherited_resolved_peak_in_band=f,full_resolved_peak_over_band=(peak/f if f>0. else None),
                             retained_band_has_no_overlap=bool(lo>=q[-1] or hi<=q[0])))
        rows.append(dict(T_anchor_GeV=S["T_visible_GeV"],source_shape_sha256=sha(p),finite_support_upper_GeV=q[-1],
                         all_resolved_scalar_peak=peak,frozen_scalar_n_GeV3=N,total_quadrature_relative=full/N-1.,
                         interval_probes=data))
    return dict(cases=rows,scope="At a future expansion probe, treat only previously resolved scalar occupation as a freely redshifted inherited component. Stable incident momentum assumed bounded by P/b. The interval contains the full possible final band but is not an asserted support equality. No new production, exterior return, redistribution or nonlinear Bose feedback is bounded.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shapes",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.shapes),indent=2,allow_nan=False))
