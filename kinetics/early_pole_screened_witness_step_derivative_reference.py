"""Check the derivative of a nonlinear screened split about a short archived witness."""
import argparse
import json

import numpy as np
from scipy.linalg import block_diag,expm,null_space,solve_triangular
from scipy.special import expit

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
from early_pole_screened_elastic_discrete_gradient_reference import approximate_jacobian,step


def run(grid,subdivision,exposures,macro):
    p=next(p for p in json.loads(PROFILES.read_text()) if p['grid']==grid)
    C=Screened(grid,1,16,max(exposures),(1,1,2),.08,1.,subdivision=subdivision)
    C,_,old,new,*_=initial_projection(p,subdivision,16,max(exposures),model=C)
    m=ElasticChannels(C,grid,1.)
    M=m.M
    Y=.03*old+.97*new
    _,F,G,*_=M.weights(Y,True)
    k=M.r[1].start
    cells=M.cells[k:]
    scale=np.repeat([c.W for c in cells],2).reshape(-1,2)
    blocks=[]
    for c,v in zip(cells,F[k:],strict=True):
        f=expit(v@c.B)
        B=(c.I*(f*(1-f)))@c.B.T
        blocks.append((B+B.T)/2)
    S=block_diag(*blocks)
    L=np.linalg.cholesky(S)
    J=S-approximate_jacobian(m,G,F,k,1.,scale)*scale.ravel()[:,None]
    cons=[]
    for group in M.r[1:]:
        z=np.zeros_like(F)
        z[group,0]=1.
        cons.append(z[k:].ravel())
    U=null_space(np.stack([*cons,G[k:].ravel()])@L)
    A=solve_triangular(L,solve_triangular(L,J,lower=True).T,lower=True).T
    B=U.T@A@U
    v=U.T@solve_triangular(L,(new-old)[k:].ravel(),lower=True)
    def theta(w):
        return solve_triangular(L.T,U@w,lower=False)
    v/=max(abs(theta(v)))
    seed=theta(v)

    def material(w,delta):
        logits=F.copy()
        logits[k:]+=delta*w.reshape(-1,2)
        out=Y.copy()
        for j,c in enumerate(cells,start=k):
            out[j]=c.I@expit(logits[j]@c.B)
        return out
    def source(y):
        z=np.r_[y.ravel(),C.ns,0.,0.]
        return float(Cold.rhs(C,z,False)[C.j0:C.j0+C.n].sum())
    def derivative(pair,eps):
        dy=(pair[1]-pair[0])[k:].ravel()/(2*eps)
        return (float(np.linalg.norm(solve_triangular(L,dy,lower=True))/np.linalg.norm(v)),
                (source(pair[1])-source(pair[0]))/(2*eps))

    out=[]
    for t in exposures:
        rows=[]
        for eps in (3e-4,1e-4):
            yl,yr=(material(seed,d*eps) for d in (-1.,1.))
            for _ in range(macro):
                yl,_=step(m,yl,t/macro)
                yr,_=step(m,yr,t/macro)
            w=expm(t*B)@v
            rows.append(dict(perturbation=eps,nonlinear_direction_norm_and_source=derivative((yl,yr),eps),
                             frozen_direction_norm_and_source=derivative((material(theta(w),-eps),material(theta(w),eps)),eps)))
        out.append(dict(exposure=t,nonlinear_macro_steps=macro,rows=rows))
    return dict(grid=grid,subdivision=subdivision,results=out,
                scope='Directional derivative of a finite conservative elastic step vs held frozen tangent; no screening refresh or slow channels.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--subdivision',type=int,default=3)
    p.add_argument('--exposures',type=float,nargs='+',default=[1e-10,1e-9,1e-8])
    p.add_argument('--macro',type=int,default=1)
    a=p.parse_args()
    if a.macro<1 or min(a.exposures)<=0.:
        p.error('Choose positive exposures and at least one macro step.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
