"""Match retained daughter number and model total hard energy at fixed local cells."""
import argparse
import gzip
import json

import numpy as np
from scipy.optimize import root
from scipy.special import expit

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_complement_retune_reference import prepare


def dressed(M,native,continuous,fixed):
    raw,bare,hints=prepare(M,native,continuous,fixed)
    daughter=np.arange(M.r[2].start,M.r[2].stop)
    donor=np.asarray(hints['complement_cells'])
    shape=M.shape
    number=lambda y:float(y[daughter,0].sum())
    bare_energy=lambda y:float(np.sum(y[daughter]*shape[daughter]))
    hard_energy=lambda y:float(M.weights(y,True)[3])
    number_ref=number(native)
    energy_ref=hard_energy(native)
    pivot=(bare_energy(native)-np.sum(raw[fixed]*shape[fixed]))/(number_ref-raw[fixed,0].sum())
    bases=[np.vstack([np.ones_like(M.cells[i].B[0]),shape[i]@M.cells[i].B-pivot])
           for i in donor]
    fd=[M.cells[i].fit(native[i])[0]@M.cells[i].B for i in donor]
    def evaluate(theta):
        y=raw.copy()
        slopes=[]
        for i,t,B in zip(donor,fd,bases,strict=True):
            cell=M.cells[i]
            f=expit(t+theta@B)
            y[i]=cell.I@f
            slopes.append((cell.I*(f*(1-f)))@B.T)
        _,_,grad,energy,_=M.weights(y,True)
        jac=np.array([np.sum([s[0] for s in slopes],axis=0),
                      np.sum([grad[i]@s for i,s in zip(donor,slopes,strict=True)],axis=0)])
        return np.array([number(y)-number_ref,energy-energy_ref]),jac,y
    theta=root(lambda t:evaluate(t)[0],np.asarray(hints['affine_logit_multipliers']),
               jac=lambda t:evaluate(t)[1],tol=1e-11).x
    defect,jac,y=evaluate(theta)
    if np.max(abs(defect))>=5e-14 or np.max(abs(theta))>=.2:
        raise ValueError('Modeled two-moment amplitude or accuracy gate failed.')
    check=1e-5
    fd_jac=np.column_stack([(evaluate(theta+v*check)[0]-evaluate(theta-v*check)[0])/(2*check)
                            for v in np.eye(2)])
    np.testing.assert_allclose(jac,fd_jac,rtol=3e-6,atol=2e-10)
    untouched=np.setdiff1d(np.arange(len(y)),donor)
    np.testing.assert_array_equal(y[untouched],raw[untouched])
    out=dict(fixed_cells=list(map(int,fixed)),donor_cells=hints['complement_cells'],
             coefficients=theta.tolist(),number_and_modeled_energy_residual=defect.tolist(),
             finite_difference_jacobian_max_gap=float(np.max(abs(jac-fd_jac))),
             bare_match_number_and_modeled_energy_residual=[number(bare)-number_ref,hard_energy(bare)-energy_ref],
             dressed_match_bare_energy_residual=bare_energy(y)-bare_energy(native),
             max_absolute_cell_moment_difference_vs_bare_match=float(np.max(abs(y-bare))))
    return raw,bare,y,out


def run(bundle):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    a,b=(branches[k] for k in ('steps','cell-conservative'))
    u,v=(q['snapshots'][0] for q in (a,b))
    n=np.asarray(u['cold_bins'])
    x,y=(np.asarray(q['hard_moments']) for q in (u,v))
    C=Cold(a['grid'],1,0,len(n),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection='steps')
    M=C.M
    np.testing.assert_array_equal(n,v['cold_bins'])
    np.testing.assert_array_equal(n,C.ns)
    np.testing.assert_array_equal(x,M.init)
    for e,q,z in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(e,q)
        np.testing.assert_array_equal(e,z)
    ed=np.asarray(M.edge_sets[2])
    ix=M.r[2].start+np.flatnonzero((ed[:-1]>=1.4138888888888888-2e-14)&(ed[1:]<=1.4666666666666666+2e-14))
    if len(ix)!=2:
        raise ValueError('Expected the original local two-cell intervention.')
    raw,bare,full,out=dressed(M,x,y,ix)
    source=lambda t:Cold.rhs(C,np.r_[t.ravel(),n,0.,0.],False)[C.j0:-2]
    ref=source(x)
    np.testing.assert_array_equal(ref,u['instantaneous_source_by_q'])
    return dict(preparation=out,initial_rate_shifts={
                    k:float(np.sum(source(z)-ref)) for k,z in
                    [('raw_local',raw),('bare_match',bare),('modeled_match',full)]},
                scope='Static same retained daughter number and modeled hard energy; cold seed fixed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    a=p.parse_args()
    print(json.dumps(run(a.bundle),indent=2,allow_nan=False))
