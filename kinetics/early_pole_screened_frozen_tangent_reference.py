"""Frozen-geometry tangent of the elastic operators, not a spectral width."""
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
from early_pole_screened_elastic_discrete_gradient_reference import approximate_jacobian


def run(grid, sector, subdivision):
    profile=next(p for p in json.loads(PROFILES.read_text()) if p['grid']==grid)
    duration=2e-8
    C=Screened(grid,1,16,duration,(1,1,2),.08,float(sector!='same'),subdivision=subdivision)
    C,archive,Y0,Y1,*_=initial_projection(profile,subdivision,16,duration,model=C)
    m=ElasticChannels(C,grid,float(sector!='unlike'))
    M=m.M
    k=M.r[1].start
    cells=M.cells[k:]
    scale=np.repeat([c.W for c in cells],2).reshape(-1,2)
    Y=.03*Y0+.97*Y1
    _,F,G,*_=M.weights(Y,True)
    case=[]
    for thermal in (False,True):
        base=F.copy()
        if thermal:
            for group,mu in zip(M.r[1:],(.01,-.02),strict=True):
                base[group]=-G[group]
                base[group,0]+=mu
        block=[]
        for c,v in zip(cells,base[k:],strict=True):
            f=expit(v@c.B)
            B=(c.I*(f*(1-f)))@c.B.T
            block.append((B+B.T)/2)
        S=block_diag(*block)
        L=np.linalg.cholesky(S)
        J=S-approximate_jacobian(m,G,base,k,1.,scale)*scale.ravel()[:,None]
        count=[]
        for group in M.r[1:]:
            z=np.zeros_like(base)
            z[group,0]=1.
            count.append(z[k:].ravel())
        energy=[]
        if sector=='same':
            for group in M.r[1:]:
                z=np.zeros_like(base)
                z[group]=G[group]
                energy.append(z[k:].ravel())
        else:
            energy.append(G[k:].ravel())
        cons=np.stack([*count,*energy])
        res=np.linalg.norm(cons@J,axis=1)/np.maximum(1e-100,np.linalg.norm(cons,axis=1)*np.linalg.norm(J))
        assert np.max(res)<2e-12
        U=null_space(cons@L)
        A=solve_triangular(L,solve_triangular(L,J,lower=True).T,lower=True).T
        B=U.T@A@U
        anorm=np.linalg.norm(B)
        eigen,vectors=np.linalg.eig(B)
        on=np.argmax(eigen.real)
        ev=vectors[:,on]
        cond=float(np.linalg.cond(S))
        spec=float(eigen[on].real)
        eigen_residual=float(np.linalg.norm(B@ev-eigen[on]*ev)/(anorm*np.linalg.norm(ev)))
        assert eigen_residual<3e-12
        maxsym=float(np.linalg.eigvalsh((B+B.T)/2)[-1])
        nonsym=np.linalg.norm(B-B.T)/max(anorm,1e-100)
        if thermal:
            assert nonsym<2e-9 and maxsym<=2e-10*anorm
        init=solve_triangular(L,(Y1-Y0)[k:].ravel(),lower=True)
        tangent=U.T@init
        lost=np.linalg.norm(init-U@tangent)/np.linalg.norm(init)
        def logits(v):
            return solve_triangular(L.T,U@v,lower=False)
        tangent/=max(abs(logits(tangent)))
        q=logits(tangent)
        def diff(v,h):
            lo,hi=base.copy(),base.copy()
            lo[k:]-=h*v.reshape(-1,2);hi[k:]+=h*v.reshape(-1,2)
            return (m.collision(G,hi,G)[0][k:].ravel()-m.collision(G,lo,G)[0][k:].ravel())/(2*h)
        target=J@q
        err=[float(np.linalg.norm(diff(q,h)-target)/np.linalg.norm(target)) for h in (3e-4,1e-4)]
        assert max(err)<5e-6
        def cold(v,h):
            def source(a):
                th=base.copy()
                th[k:]+=a*v.reshape(-1,2)
                Yout=Y.copy()
                for j,c in enumerate(cells,start=k):
                    Yout[j]=c.I@expit(th[j]@c.B)
                z=np.r_[Yout.ravel(),C.ns,0.,0.]
                return Cold.rhs(C,z,False)[C.j0:C.j0+C.n].sum()
            return float((source(h)-source(-h))/(2*h))
        cold0=cold(q,3e-4)
        abs_readout_error=abs(cold(q,1e-4)-cold0)
        rows=[]
        for time in (0.,1e-10,3e-10,1e-9,3e-9,1e-8,3e-8,1e-7):
            u=expm(time*B)@tangent
            rows.append(dict(exposure=time,weighted_moment_norm_ratio=float(np.linalg.norm(u)/np.linalg.norm(tangent)),
                             overlap_ratio=float(u@tangent/(tangent@tangent)),
                             scalar_source_directional_derivative=cold(logits(u),3e-4)))
        case.append(dict(common_affine_control=thermal,normalized_conservation_defects=res.tolist(),
                         jacobian_directional_errors=err,relative_constraint_projection=lost,
                         relative_nonsymmetric_norm=float(nonsym),maximum_symmetric_eigenvalue=maxsym,
                         susceptibility_condition_number=cond,spectral_abscissa=spec,
                         leading_eigenpair_relative_residual=eigen_residual,
                         initial_readout_absolute_gate=abs_readout_error,
                         initial_cold_source_derivative=cold0,rows=rows))
    return dict(grid=grid,sector=sector,subdivision=subdivision,cases=case,
                scope='Only elasticity and its frozen graph, occupations and energy metric; not spectral damping or a live nonlinear trajectory.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('grid',choices=('fine','third'))
    ap.add_argument('--sector',choices=('same','unlike','both'),default='both')
    ap.add_argument('--subdivision',type=int,default=3)
    args=ap.parse_args()
    print(json.dumps(run(**vars(args)),indent=2,allow_nan=False))
