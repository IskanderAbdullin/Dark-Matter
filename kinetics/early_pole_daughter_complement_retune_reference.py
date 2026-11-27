"""Retune only the unswapped daughter cells to match native number and bare energy."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import root
from scipy.special import expit

from early_pole_archived_shared_cold_window_reference import Cold


def prepare(M,native,continuous,fixed_ids):
    full=np.arange(M.r[2].start,M.r[2].stop)
    fixed=np.asarray(fixed_ids)
    complement=np.setdiff1d(full,fixed,assume_unique=True)
    if not 0<len(complement)<len(full):
        raise ValueError('Require a nonempty proper daughter subset.')
    swapped=native.copy()
    swapped[fixed]=continuous[fixed]
    weights=M.shape
    moments=lambda y,ids:np.array([y[ids,0].sum(),np.sum(y[ids]*weights[ids])])
    target=moments(native,full)-moments(swapped,fixed)
    pivot=target[1]/target[0]
    bases=[np.vstack([np.ones_like(M.cells[i].B[0]),weights[i]@M.cells[i].B-pivot])
           for i in complement]
    theta=[M.cells[i].fit(native[i])[0]@M.cells[i].B for i in complement]
    def eval_at(x):
        trial=swapped.copy()
        J=np.zeros((2,2))
        for i,B,t in zip(complement,bases,theta,strict=True):
            cell=M.cells[i]
            f=expit(t+x@B)
            trial[i]=cell.I@f
            diff=(cell.I*(f*(1-f)))@B.T
            J+=np.vstack([diff[0],weights[i]@diff])
        return moments(trial,complement)-target,J,trial
    fit=root(lambda x:eval_at(x)[0],np.zeros(2),jac=lambda x:eval_at(x)[1],tol=1e-11)
    _,jac,trial=eval_at(fit.x)
    defect=moments(trial,full)-moments(native,full)
    if np.max(abs(defect))>=5e-14 or np.max(abs(fit.x))>=.2:
        raise ValueError('Retune failed the fixed-moment or local-amplitude gate: '
                         f'defect={defect.tolist()}, multipliers={fit.x.tolist()}')
    np.testing.assert_array_equal(trial[fixed],continuous[fixed])
    return swapped,trial,dict(fixed_cells=fixed.tolist(),complement_cells=complement.tolist(),
                              number_and_bare_energy_defect=defect.tolist(),
                              original_swap_moment_defect=(moments(swapped,full)-moments(native,full)).tolist(),
                              affine_logit_multipliers=fit.x.tolist(),
                              jacobian_condition=float(np.linalg.cond(jac)))


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    a,b=(branches[k] for k in ('steps','cell-conservative'))
    Y0,Y1=(np.asarray(x['snapshots'][0]['hard_moments']) for x in (a,b))
    N=np.asarray(a['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(N,b['snapshots'][0]['cold_bins'])
    C=Cold(a['grid'],1,0,len(N),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    M=C.M
    np.testing.assert_array_equal(M.init,Y0)
    for e,x,y in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(e,x)
        np.testing.assert_array_equal(e,y)
    d=M.r[2]
    edge=np.asarray(M.edge_sets[2])
    q=np.flatnonzero((edge[:-1]>=1.4138888888888888-2e-14)&(edge[1:]<=1.4666666666666666+2e-14))
    if len(q)!=2:
        raise ValueError('Expected two original cold-active daughter cells.')
    mode=lambda y:Cold.rhs(C,np.r_[y.ravel(),N,0.,0.],False)[C.j0:-2]
    out={}
    ref=mode(Y0)
    for key,ids in (('local',d.start+q),('outside',np.setdiff1d(np.arange(d.start,d.stop),d.start+q))):
        try:
            p,r,check=prepare(M,Y0,Y1,ids)
        except ValueError as exc:
            out[key]={'rejected':str(exc)}
            continue
        old,new=(mode(x) for x in (p,r))
        out[key]=dict(check,original_shift=float(old.sum()-ref.sum()),
                      retuned_shift=float(new.sum()-ref.sum()),
                      retuned_mode_shift=(new-ref).tolist())
    return dict(native_rate=float(ref.sum()),cases=out,
                scope='Fixed bare daughter energy and number with complementary controlled FD retuning; initial source only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
