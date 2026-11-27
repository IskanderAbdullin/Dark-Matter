"""Conditional hard-fermion mixing at exactly fixed energy and species counts."""
import argparse
import json

import numpy as np
from scipy.special import expit
from early_pole_archived_shared_cold_window_reference import Cold


class Mixing(Cold):
    def __init__(self,grid,wide,n,rate,mode,exposure):
        super().__init__(grid,wide,0,n,True,exposure=exposure)
        self.rate=rate
        self.mode=mode
        self.beg=self.M.r[1].start
        self.split=self.M.r[2].start-self.beg
        k=len(self.M.cells)-self.beg
        self.c=np.zeros((k,2,2))
        self.c[:self.split,0,0]=1.
        self.c[self.split:,0,1]=1.

    def mixing(self,y):
        Y=y[:self.j0].reshape(-1,2)
        _,F,E,_,_=self.M.weights(Y,True)
        F,E=F[self.beg:],E[self.beg:]
        metric=[]
        for C,Q in zip(self.M.cells[self.beg:],F,strict=True):
            O=expit(Q@C.B)
            metric.append((C.I*(O*(1-O)))@C.B.T)
        metric=np.asarray(metric)
        if self.mode=='joint':
            constraint=np.concatenate([self.c,E[...,None]],axis=2)
        elif self.mode=='separate':
            e0,e1=np.zeros_like(E),np.zeros_like(E)
            e0[:self.split]=E[:self.split]
            e1[self.split:]=E[self.split:]
            constraint=np.concatenate([self.c,e0[...,None],e1[...,None]],axis=2)
        else:
            raise ValueError(self.mode)
        n=F.shape[0]*2
        C=constraint.reshape(n,-1)
        M=lambda a:np.einsum('ijk,ik...->ij...',metric,a.reshape((len(metric),2)+a.shape[1:])).reshape(a.shape)
        S=-F.ravel()
        K=C.T@M(C)
        G=np.sqrt(np.diag(K))
        if min(G)<=0.:raise ValueError('degenerate thermodynamic constraint')
        Q=C/G
        driver=S-Q@np.linalg.solve(Q.T@M(Q),Q.T@M(S))
        D=self.rate*M(driver)
        close=C.T@D
        ent=S@D
        pos=self.rate*driver@M(driver)
        assert max(abs(close))<2e-10*max(1.,np.linalg.norm(D))
        assert abs(ent-pos)<1e-9*max(1.,abs(ent)) and pos>=0.
        return D,ent,close

    def rhs(self,y,pair):
        f=super().rhs(y,pair)
        if self.rate:
            R,H,C=self.mixing(y)
            f[2*self.beg:self.j0]+=R
            f[-2]+=H
            f[-1]+=sum(C[2:])
            # At fixed y the drift is separately audited; actual finite
            # energy conservation is checked in the enclosing integration.
        return f


def static_witness(M):
    ini=np.r_[M.M.init.ravel(),M.ns,0.,0.]
    d,s,c=M.mixing(ini)
    E=M.M.weights(M.M.init,True)[2][M.beg:]
    D=d.reshape(-1,2)
    heat=[np.sum(E[:M.split]*D[:M.split]),np.sum(E[M.split:]*D[M.split:])]
    residual=[]
    for common in (True,False):
        Y=M.M.init.copy()
        for k in range(48):
            old=Y.copy()
            E=M.M.weights(Y,True)[2]
            for j,C in enumerate(M.M.cells[M.beg:],M.beg):
                b=1. if j<M.M.r[2].start or common else 1.6
                mu=-.3 if j<M.M.r[2].start else -.65
                Y[j]=C.I@expit(-b*(E[j,0]+E[j,1]*C.u-mu))
            if np.max(abs(Y-old))<3e-15:
                break
        z=np.r_[Y.ravel(),M.ns,0.,0.]
        D,S,cons=M.mixing(z)
        mx=max(abs(D))
        if common or M.mode=='separate': assert mx<3e-10
        else: assert S>1e-5
        residual.append(dict(common_temperature=common,response_inf=mx,entropy_rate=S))
    return dict(initial_energy_exchange_per_reference_log_a=heat,initial_entropy_production=s,
                initial_constraint_defect=c.tolist(),manufactured_thermal=residual)


def run(grid,wide,n,rate,mode,exposure):
    M=Mixing(grid,wide,n,rate,mode,exposure)
    o=M.run()
    o.update(mixing_rate=rate,mixing_mode=mode,static_audit=static_witness(M),scope='Conditional entropy-gradient mixing, '
             'not a calibrated physical elastic operator. Photons not mixed in this added term.')
    return o


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',default='fine',choices=['fine','third'])
    p.add_argument('--wide',default=1,type=int,choices=[0,1,2])
    p.add_argument('--n',default=8,type=int)
    p.add_argument('--rate',default=500.,type=float)
    p.add_argument('--mode',choices=['joint','separate'],default='joint')
    p.add_argument('--exposure',type=float,default=.0003)
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.wide,a.n,a.rate,a.mode,a.exposure),indent=2))
