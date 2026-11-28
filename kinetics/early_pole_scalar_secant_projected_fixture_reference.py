"""A manufactured positive Hartree-secant projection, not a new physical rate."""
import json

import numpy as np
from scipy.optimize import root
from scipy.sparse import coo_matrix
from scipy.special import roots_legendre

from early_pole_scalar_elastic_stencil_reference import Elastic
from early_pole_scalar_hartree_discrete_gradient_reference import energy, entropy, gap, gradient


def project(graph,secant):
    g=secant-secant[0]
    i,j,l1,h1,l2,h2=graph.j.T
    a,b=graph.v[:,3],graph.v[:,5]
    d1,d2=g[h1]-g[l1],g[h2]-g[l2]
    target=g[i]+g[j]-g[l1]-g[l2]
    delta=target-a*d1-b*d2
    up=delta>=0.
    u1=np.where(up,1-a,-a)
    u2=np.where(up,1-b,-b)
    slack=d1*u1+d2*u2
    capable=(target>=-2e-15*np.max(g))&(target<=d1+d2+2e-15*np.max(g))
    factor=np.divide(delta,slack,out=np.zeros_like(delta),where=(slack!=0.)&capable)
    factor=np.where(capable,np.clip(factor,0.,1.),0.)
    aa,bb=np.clip(a+factor*u1,0.,1.),np.clip(b+factor*u2,0.,1.)
    v=np.c_[-np.ones((len(i),2)),1-aa,aa,1-bb,bb]
    v[~capable]=0.
    rows=np.repeat(np.arange(len(i)),6)
    S=coo_matrix((v.ravel(),(rows,graph.j.ravel())),shape=(len(i),len(g))).tocsr()
    geom=np.maximum(abs((aa-a)*(graph.p[h1]-graph.p[l1])),
                    abs((bb-b)*(graph.p[h2]-graph.p[l2])))
    return S,capable,geom


def case(edge,amplitude):
    q=np.cbrt((edge[1:]**3+edge[:-1]**3)/2)
    V=np.diff(edge**3)/(6*np.pi**2)
    lam,ub=7.65175e-8,(1e-4)**2
    n=np.exp(amplitude*np.sin(3*np.linspace(0,1,len(q))))*V**.4
    n*=.017/n.sum()
    u=gap(q,n,ub,lam)
    K=Elastic(q,np.sqrt(u),edge,8,lam)
    x,w=roots_legendre(20)
    R,A=K.reaction(n/V)
    gain=np.divide(R,A,out=np.zeros_like(R),where=A!=0.)
    assert np.min(gain)>=0.
    norm=1/max(abs((K.S.T@R)/n))
    h=.08*norm
    pivot=int(np.argmax(n))
    free=np.flatnonzero(np.arange(len(n))!=pivot)
    total=float(np.sum(n))
    def evaluate(z):
        t=n.copy()
        t[free]*=np.exp(z)
        t[pivot]=total-np.sum(t[free])
        if t[pivot]<=0.:
            raise ValueError("A trial left the positive fixed-number simplex.")
        v=gap(q,t,ub,lam)
        E=gradient(q,n,t,np.sqrt(u),np.sqrt(v),lam)
        S,keep,geo=project(K,E)
        p=n[:,None]+.5*(t-n)[:,None]*(1+x)
        ds=.5*np.log1p(V[:,None]/p)@w
        F=S@ds
        change=h*(S.T@(gain*F))
        return t,v,E,S,keep,geo,F,change

    def residual(z):
        t,*_,change=evaluate(z)
        return ((t-n-change)/n)[free]

    m=root(residual,np.zeros(len(free)),tol=4e-13)
    t,v,E,S,keep,geo,F,change=evaluate(m.x)
    full_residual=float(np.max(abs((t-n-change)/n)))
    assert full_residual<3e-12
    de=energy(q,t,v,ub,lam)-energy(q,n,u,ub,lam)
    dS=float(np.sum(entropy(t,V)-entropy(n,V)))
    entropy_ledger=h*np.dot(gain,F**2)
    es=float(np.max(abs(S@(E-E[0]))))
    assert abs(de)<5e-16 and abs(t.sum()-n.sum())<2e-14
    assert dS>0 and abs(dS-entropy_ledger)<1e-9*abs(dS)
    assert es<5e-15*np.max(E-E[0])
    return dict(q_support=[edge[0],edge[-1]],cells=len(q),amplitude=amplitude,
                relative_number_defect=float(abs(t.sum()-n.sum())/n.sum()),
                finite_energy_defect=float(de),entropy_change=dS,entropy_ledger=float(entropy_ledger),
                max_projection_secant_defect=es,projected_events=int(keep.sum()),
                rejected_event_fraction=float(np.mean(~keep)),
                max_per_leg_momentum_magnitude_change_over_qmax=float(max(geo)/max(q)),
                max_scaled_nonlinear_residual=full_residual)


if __name__=='__main__':
    print(json.dumps([case(np.geomspace(a,b,17),s)
                      for a,b in ((1e-6,1e-5),(1e-6,3e-3)) for s in (.7,1.6)],indent=2,allow_nan=False))
