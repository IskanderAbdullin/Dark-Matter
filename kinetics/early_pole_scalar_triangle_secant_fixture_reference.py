"""Manufactured Hartree-secant map restricted to feasible momentum triangles."""
import json

import numpy as np
from scipy.optimize import root
from scipy.sparse import coo_matrix
from scipy.special import roots_legendre

from early_pole_scalar_elastic_stencil_reference import Elastic
from early_pole_scalar_hartree_discrete_gradient_reference import energy, entropy, gap, gradient


def project(K,g,panel):
    n=len(g)
    i,j,L1,H1,L2,H2=K.j.T
    w1,w2=K.v[:,3],K.v[:,5]
    d=np.arange(-panel,2+panel)
    A=np.clip(L1[:,None]+d,0,n-1)
    B=np.clip(L2[:,None]+d,0,n-1)
    A=np.broadcast_to(A[:,:,None],(len(i),len(d),len(d))).reshape(len(i),-1)
    B=np.broadcast_to(B[:,None,:],(len(i),len(d),len(d))).reshape(len(i),-1)
    u,v=K.p[A],K.p[B]
    P=K.lab[:,None]
    eps=4e-15*max(K.p)
    allowed=(np.abs(u-v)<=P+eps)&(u+v>=P-eps)
    E=g-g[0]
    e=E[A]+E[B]
    assert np.all(np.diff(E)>0.)
    target=E[i]+E[j]
    mn=np.where(allowed,e,np.inf).min(axis=1)
    mx=np.where(allowed,e,-np.inf).max(axis=1)
    take=(target>=mn-1e-15*max(E))&(target<=mx+1e-15*max(E))
    prior=(np.where(A==L1[:,None],1-w1[:,None],np.where(A==H1[:,None],w1[:,None],0.))
           *np.where(B==L2[:,None],1-w2[:,None],np.where(B==H2[:,None],w2[:,None],0.)))
    prior=np.where(allowed,prior,0.)
    min_idx=np.where(allowed,e,np.inf).argmin(axis=1)
    max_idx=np.where(allowed,e,-np.inf).argmax(axis=1)
    idx=np.arange(len(i))
    use=prior.sum(axis=1)
    r=np.divide(prior,use[:,None],where=use[:,None]>0,out=np.zeros_like(prior))
    r[idx[use==0.],min_idx[use==0.]]=1.
    mean=(r*e).sum(axis=1)
    choose=np.where(target>=mean,max_idx,min_idx)
    cap=e[idx,choose]-mean
    weight=np.clip(np.divide(target-mean,cap,where=take&(cap!=0.),out=np.zeros(len(i))),0.,1.)
    r*=1-weight[:,None]
    r[idx,choose]+=weight
    r[~take]=0.
    positions=np.c_[i,j,A,B]
    values=np.c_[np.where(take,-1.,0.),np.where(take,-1.,0.),r,r]
    S=coo_matrix((values.ravel(),(np.repeat(idx,values.shape[1]),positions.ravel())),
                 shape=(len(i),len(g))).tocsr()
    support=(r>0.)
    geo=np.maximum(0.,np.maximum(abs(u-v)-P,P-u-v))
    return S,take,float(np.max(geo[support],initial=0.))/max(K.p)


def one(edge,panel):
    q=np.cbrt((edge[1:]**3+edge[:-1]**3)/2)
    V=np.diff(edge**3)/(6*np.pi**2)
    lam,ub=7.65175e-8,1e-8
    n=V**.4*np.exp(1.25*np.sin(3*np.linspace(0,1,len(q))))
    n*=.017/n.sum()
    x,w=roots_legendre(20)
    u=gap(q,n,ub,lam)
    K=Elastic(q,np.sqrt(u),edge,8,lam,keep_lab=True)
    R,A=K.reaction(n/V)
    M=np.divide(R,A,out=np.zeros_like(R),where=A!=0.)
    assert np.min(M)>=0.
    tau=.05/np.max(abs((K.S.T@R)/n))
    pivot=np.argmax(n)
    free=np.flatnonzero(np.arange(len(n))!=pivot)
    count=n.sum()
    def evaluate(z):
        nn=n.copy()
        nn[free]*=np.exp(z)
        nn[pivot]=count-nn[free].sum()
        if nn[pivot]<=0.:
            raise ValueError("nonpositive pivot")
        v=gap(q,nn,ub,lam)
        g=gradient(q,n,nn,np.sqrt(u),np.sqrt(v),lam)
        d=.5*np.log1p(V[:,None]/(n[:,None]+.5*(nn-n)[:,None]*(1+x)))@w
        S,keep,geom=project(K,g,panel)
        F=S@d
        return nn,v,g,S,keep,geom,F,tau*(S.T@(M*F))
    def eq(z):
        nn,*_,flux=evaluate(z)
        return ((nn-n-flux)/n)[free]
    solve=root(eq,np.zeros(len(free)),tol=4e-13)
    nn,v,g,S,keep,geom,F,flux=evaluate(solve.x)
    res=float(np.max(abs((nn-n-flux)/n)))
    de=energy(q,nn,v,ub,lam)-energy(q,n,u,ub,lam)
    ds=float(np.sum(entropy(nn,V)-entropy(n,V)))
    sl=tau*np.dot(M,F*F)
    error=float(np.max(abs(S@(g-g[0]))))
    assert ds>0. and abs(ds-sl)<1e-9*abs(ds)
    assert error<8e-15*max(g-g[0]) and geom<5e-15
    assert abs(de)<4e-16 and res<3e-12
    return dict(q_range=[edge[0],edge[-1]],extra_outgoing_neighbors=panel,accepted=int(keep.sum()),
                geometrical_input_count=len(keep),retained_born_weight=float(K.weight[keep].sum()/K.weight.sum()),
                triangle_violation_over_qmax=geom,scalar_number_error=float(nn.sum()-count),
                finite_energy_error=de,entropy_change=ds,independent_entropy_ledger=float(sl),
                maximum_secant_direction_defect=error,scaled_nonlinear_residual=res)


if __name__=='__main__':
    print(json.dumps([one(np.geomspace(a,b,17),panel)
                      for a,b in ((1e-6,1e-5),(1e-6,3e-3)) for panel in (0,1)],indent=2,allow_nan=False))
