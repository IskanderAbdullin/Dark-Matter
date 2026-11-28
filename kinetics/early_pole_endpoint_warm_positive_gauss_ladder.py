"""Compare positive q²-volume Gauss orders on the warm frozen collision source."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_source_cells import local


def orthogonal_volume_gauss(a,b,p):
    """Use bounded orthogonal polynomials instead of nearly collinear monomials."""
    z,W=roots_legendre(p+2)
    t=(1+z)/2
    P=np.polynomial.legendre.legvander(z,p-1)*np.sqrt(2*np.arange(p)+1)
    A,D=a/b,(b-a)/b
    q=[]
    g=[]
    for u,v in zip(A,D,strict=True):
        rho=W*(u+v*t)**2
        rho/=rho.sum()
        B=(P.T*rho)@P
        X=(P.T*(rho*t))@P
        nodes,vec=eigh(X,B)
        weight=(B[0]@vec)**2
        q.append(nodes)
        g.append(weight)
    return np.array(q),np.array(g)


def volume_gauss(a,b,p):
    """Work on u=(q-a)/(b-a): H and X are the volume Gram and multiplication matrices."""
    a,b=np.asarray(a),np.asarray(b)
    if a.ndim!=1 or a.shape!=b.shape or not np.all(np.isfinite([a,b])) or np.min(a)<0 or np.any(a>=b) or p not in range(2,13):
        raise ValueError('Use positive finite shells and Gauss orders 2 through 12')
    A,H=a/b,(b-a)/b
    t=np.arange(2*p,dtype=float)
    I=(A[:,None]**2/(1+t)+(2*A*H)[:,None]/(2+t)+H[:,None]**2/(3+t))
    I/=I[:,:1].copy()
    if p>=7:
        u,w=orthogonal_volume_gauss(a,b,p)
    else:
        # Preserve the exact lower-order archives and the originally checked rule.
        ix=np.arange(p)[None,:]+np.arange(p)[:,None]
        U=[]
        W=[]
        for z in I:
            B=z[ix]
            u,v=eigh(z[ix+1],B)
            w=(B[0]@v)**2
            U.append(u)
            W.append(w)
        u,w=np.array(U),np.array(W)
    if np.any(u<=0) or np.any(u>=1) or np.any(w<=0):
        raise ArithmeticError('A node or weight is not interior and positive')
    err=np.max(abs(np.einsum('ij,ijk->ik',w,u[:,:,None]**t)/I-1))
    if err>2e-12:
        raise ArithmeticError('The declared volume moments did not close')
    vol=(b-a)*(a*a+a*b+b*b)/(6*np.pi**2)
    return a[:,None]+(b-a)[:,None]*u,w*vol[:,None],err


class WarmGauss(ScalarSubshellChart):
    def __init__(self,bundle,p):
        super().__init__(bundle,SPECS['fine'],2,4,1e-5)
        a,b=self.ds.ed[:-1],self.ds.ed[1:]
        take=a>=.01
        self.old_count=self.ns
        q,w,err=volume_gauss(a,b,p)
        self.polynomial_error=err
        w*=self.ds.w[:,None]/w.sum(axis=1)[:,None]
        self.repeats=np.where(take,p,1)
        sel=np.arange(p)[None,:]<self.repeats[:,None]
        q[~take,0]=self.ds.p[~take]
        w[~take,0]=self.ds.w[~take]
        self.ds.p=q[sel]
        self.ds.w=w[sel]
        self.ds.n=self.ns=len(self.ds.p)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.initial=self.transfer_step(self.initial)

    def transfer_step(self,y):
        if y.shape!=(self.old_count+self.np+self.nf,):
            raise ArithmeticError('One scalar occupation per original cell required')
        return np.r_[np.repeat(y[:self.old_count],self.repeats),y[self.old_count:]]


def run(card,ref):
    x=card['x_end']
    a,y=unpack(card,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    warm=ed[:-1]>=.01
    if (ref['x'],ref['mass_held'],ref['clock_held'])!=(x,m,H) or ref['warm_native_edges']!=ed[np.searchsorted(ed,.01):].tolist():
        raise ArithmeticError('Paired reference is not the same support, state, mass and clock')
    V={r['rule']:np.array(r['warm_cell_forward_inverse_net']) for r in ref['cases']}
    if not ref['mesh_reference']['global_L1_over_total']:
        raise ArithmeticError('The reference must pass the aggregate L1 resolution gate')
    total=fsum(V['Volume64'][:,2])
    if total<=0 or abs(total/ref['ref_net_pressure']-1)>2e-14:
        raise ArithmeticError('Reference changed normalization')
    out=[]
    for p in (2,3,4):
        d=WarmGauss(load('fine')[0],p)
        u=d.transfer_step(y)
        q=d.ds.p
        P=np.where(q>=ref['warm_floor'],q*q/np.hypot(q,m),0.)
        if any(not np.array_equal(getattr(a,k),getattr(d,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Charge chart changed')
        count=np.add.reduceat(d.ds.w*u[:d.ns],np.r_[0,np.cumsum(d.repeats)[:-1]])
        n=a.ds.w*y[:a.ns]
        if np.max(abs(count-n)/np.maximum(n,1e-105))>3e-14:
            raise ArithmeticError('Native count has changed on lift')
        l,h=[local(d,u,P,x,m,H,ed,k)[warm] for k in (4,8)]
        err=np.max(abs(l-h)/np.maximum(h[:,0,None]+h[:,1,None],1e-100))
        if err>3e-6:
            raise ArithmeticError('Radial order not resolved')
        if p==2 and np.max(abs(V['Gauss2']-h)/np.maximum(h[:,0,None]+h[:,1,None],1e-100))>3e-12:
            raise ArithmeticError('General moment rule did not reproduce the analytic two-point case')
        j=h[:,2]-V['Volume64'][:,2]
        out.append(dict(moment_order=p,nodes=d.ns,relative_polynomial_error=d.polynomial_error,
                        radial_relative_error=err,relative_signed=fsum(j)/total,
                        relative_unsigned=fsum(abs(j))/total))
    return dict(x=x,orders=out,scope='Frozen warm affine moments on a finite retained chart. Volume64 is a compared reference, not a universal cell-wise error bound.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('S4',type=Path)
    a.add_argument('reference',type=Path)
    s=a.parse_args()
    print(json.dumps(run(read(s.S4),read(s.reference)),indent=2,allow_nan=False))
