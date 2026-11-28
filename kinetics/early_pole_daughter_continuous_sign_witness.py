"""Continuous piecewise-linear daughter witnesses at fixed modeled moments."""
import argparse
import json

import numpy as np
from scipy.optimize import linprog
from scipy import sparse
from scipy.special import expit

from early_pole_daughter_frozen_moment_envelope import cell_objects
from early_pole_daughter_local_fraction_curvature_attribution import initial


def run(archive,sub):
    C,n,(y0,_,y1),_,depth=initial(archive)
    M=C.M
    e,f=M.weights(y0,True)[:2]
    pole=C.scalar(n)
    active=[j for j in range(M.r[2].start,M.r[2].stop)
            if M.cells[j].lo>=1.4138888888888888-2e-14 and M.cells[j].hi<=1.4666666666666666+2e-14]
    data=[]
    start=0
    right=None
    for j in active:
        B=M.cells[j]
        sample,*_=cell_objects(C,j,e,f,sub,pole,n,True)
        lo,hi=sample.radial_bounds
        if (right is not None and abs(lo[0]-right)>5e-13) or np.max(abs(lo[1:]-hi[:-1]))>5e-13:
            raise ValueError('Separated subcell intervals would invalidate the continuous witness.')
        right=hi[-1]
        start_node=start
        start+=len(lo)
        w,func,K,p=sample(24)
        xi=(p-lo[:,None])/(hi-lo)[:,None]
        if np.any((xi<0.)|(xi>1.)):
            raise ValueError('Interpolation weights are not positive.')
        def to_nodes(a):
            l=np.sum(w*a*(1.-xi),axis=-1)
            r=np.sum(w*a*xi,axis=-1)
            return np.r_[l,0.]+np.r_[0.,r]
        A=np.vstack([to_nodes(v) for v in func])
        Y=[]
        for y in (y0,y1):
            z=B.fit(y[j])[0]
            F=expit(z[0]+z[1]*B.basis(p))
            Y.append(np.r_[y[j],np.sum(w*func[2]*F)])
        data.append((start_node,start+1,A,Y,to_nodes(K)))
    nvar=start+1
    K=np.zeros(nvar)
    A=np.zeros((3*len(data),nvar))
    Y=np.zeros((2,len(A)))
    for i,(l,h,aa,m,k) in enumerate(data):
        norm=aa[0].sum()
        A[3*i:3*i+3,l:h]=aa/norm
        Y[:,3*i:3*i+3]=np.asarray(m)/norm
        K[l:h]+=k
    G=A.T@np.linalg.inv(A@A.T)
    vmax=max(K.sum(),1e-7)
    ends=[]
    for m in Y:
        ones=sparse.csr_matrix(np.ones((nvar,1)))
        mat=sparse.vstack((sparse.hstack((sparse.eye(nvar),ones)),
                           sparse.hstack((-sparse.eye(nvar),ones))),format='csc')
        a=linprog(np.r_[np.zeros(nvar),-1.],A_ub=mat,b_ub=np.r_[np.ones(nvar),np.zeros(nvar)],
                  A_eq=np.c_[A,np.zeros(len(A))],b_eq=m,bounds=(0.,1.),method='highs')
        if not a.success:
            raise ValueError('No interior anchor: '+a.message)
        anchor=a.x[:-1]
        anchor=anchor+G@(m-A@anchor)
        if not np.all((anchor>0.)&(anchor<1.)):
            raise ValueError('No positive interior nodal anchor.')
        pair=[]
        for sgn in (1.,-1.):
            r=linprog(sgn*K/vmax,A_eq=A,b_eq=m,bounds=(0.,1.),method='highs')
            if not r.success:
                raise ValueError(r.message)
            for sk in (1e-7,1e-6,1e-5,1e-4,.001,.01):
                y=(1.-sk)*r.x+sk*anchor
                y+=G@(m-A@y)
                if np.all((y>=0.)&(y<=1.)):
                    break
            else:
                raise ValueError('Pauli-positive corrected continuous witness missing.')
            defect=float(max(abs(A@y-m)))
            if defect>5e-13:
                raise ValueError('Nonmatching exact occupation moments.')
            pair.append(dict(influence=float(K@y),normalized_moment_residual=defect,
                             lower_occupation=float(min(y)),upper_occupation=float(max(y))))
        ends.append(pair)
    n0,m0=ends
    return dict(scope='Continuous only over selected aperture, not derivative-smooth and not evolved.',
                hard_depth=depth,cold_modes=len(n),subdivision=sub,independent_nodes_per_arm=nvar,
                checked_paired_witness_range=[n0[0]['influence']-m0[1]['influence'],n0[1]['influence']-m0[0]['influence']],
                influence_certificates=ends)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=12)
    a=p.parse_args()
    print(json.dumps([run(f,a.subdivision) for f in a.archives],indent=2,allow_nan=False))
