"""Measure archived hard-reconstruction contrasts by swapping one charged input at a time."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import root
from scipy.special import expit

from early_pole_archived_shared_cold_window_reference import Cold


def retune(M,y,reference,span):
    """Match total retained number and bare energy using two FD logit multipliers."""
    ids=list(range(span.start,span.stop))
    e=M.shape[span]
    points=[e0@cell.B for e0,cell in zip(e,(M.cells[i] for i in ids),strict=True)]
    fits=[M.cells[i].fit(y[i])[0] for i in ids]
    target=np.array([reference[span,0].sum(),np.sum(e*reference[span])])
    pivot=float(target[1]/target[0])

    def state(x):
        fields=[]
        jac=np.zeros((2,2))
        for i,mu,theta in zip(ids,points,fits,strict=True):
            c=M.cells[i]
            B=np.vstack([np.ones_like(mu),mu-pivot])
            f=expit(theta@c.B+x@B)
            fields.append(c.I@f)
            deriv=(c.I*(f*(1-f)))@B.T
            jac+=np.r_[deriv[0:1],(M.shape[i]@deriv)[None,:]]
        V=np.asarray(fields)
        value=np.array([V[:,0].sum(),np.sum(V*e)])
        return value-target,jac,V
    fit=root(lambda x:state(x)[0],np.zeros(2),jac=lambda x:state(x)[1],tol=1e-11)
    dr,_,V=state(fit.x)
    assert max(map(abs,dr))<5e-14 and np.max(abs(fit.x))<.2
    new=y.copy()
    new[span]=V
    return new,dict(target_number_and_bare_energy=target.tolist(),defect=dr.tolist(),
                    affine_logit_multipliers=fit.x.tolist(),pivot=pivot)


def read(bundle):
    with gzip.open(bundle,'rt') as f:
        return json.load(f)['branches']


def run(bundle):
    raw=read(bundle)
    a,b=(raw[p] for p in ('steps','cell-conservative'))
    sa,sb=(x['snapshots'][0] for x in (a,b))
    N=np.asarray(sa['cold_bins'])
    np.testing.assert_array_equal(N,sb['cold_bins'])
    for x,z in zip(a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,z)
    C=Cold(a['grid'],1,0,len(N),True,exposure=a['duration'],reference_aperture_grid=a['aperture_grid'],
           archive_projection='steps')
    for x,y in zip(C.M.edge_sets,a['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    x,z=(np.asarray(s['hard_moments']) for s in (sa,sb))
    M=C.M
    np.testing.assert_array_equal(x,M.init)
    np.testing.assert_array_equal(x[M.r[0]],z[M.r[0]])
    p,d=M.r[1:]
    P=x.copy();P[p]=z[p]
    D=x.copy();D[d]=z[d]
    Q,match=retune(M,P,x,p)
    R=Q.copy();R[d]=z[d]
    cases={}
    for label,S in [('native',x),('parent_only',P),('daughter_only',D),('both',z),
                    ('matched_parent_only',Q),('matched_parent_with_daughter',R)]:
        out=Cold.rhs(C,np.r_[S.ravel(),N,0.,0.],False)
        r=out[C.j0:-2]
        _,_,_,hard,s=M.weights(S,True)
        cases[label]=dict(total_cold_rate=float(r.sum()),cold_rate_each_mode=r.tolist(),
                          retained_hard_counts=[float(S[k,0].sum()) for k in M.r],
                          retained_bare_parent_and_daughter_energy=[float(np.sum(M.shape[k]*S[k])) for k in (p,d)],
                          hard_model_energy=float(hard),hard_entropy=float(s),exchange_work=float(out[-1]))
    np.testing.assert_allclose(cases['native']['cold_rate_each_mode'],sa['instantaneous_source_by_q'],rtol=1e-12,atol=1e-14)
    err=np.asarray(cases['both']['cold_rate_each_mode'])-sb['instantaneous_source_by_q']
    v={k:x['total_cold_rate'] for k,x in cases.items()}
    def interact(s):
        return v[s[1]]-v[s[0]]-v['daughter_only']+v['native']
    return dict(grid=a['grid'],scope='Held common supports and common native exterior; prepared-state instantaneous rates only.',
                swap_cases=cases,continuous_replay_absolute_max_mode_gap=float(np.max(abs(err))),
                parent_retuning=match,interaction_checks=dict(raw=interact(('parent_only','both')),
                matched=interact(('matched_parent_only','matched_parent_with_daughter'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
