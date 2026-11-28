"""Joint decay and hard pair channels sharing one conservative moment state."""
import argparse
import json
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

from early_pole_live_pauli_gauge_moment_decay_reference import Tool
from early_pole_live_reciprocal_pair_moments_reference import Pair


class Joint(Pair):
    def __init__(self,g,n,scale,edges=None):
        super().__init__(g,n,edges)
        self.pair_scale=scale
        self.start=self.r[1].start
        self.scalarq=4.1e-6
        self.Ns=.012
        self.Ws=.0002
        m=.00086724
        self.lam=2*m*m*np.hypot(m,self.scalarq)/self.Ns
        self.decay=SimpleNamespace(n=self.r[2].start-self.start,
                                   cells=self.cells[self.start:],
                                   q=self.scalarq,nod=self.z,w=self.w)

    def scalar(self,n):
        q=self.scalarq
        E=brentq(lambda e:e*(e*e-q*q)-.5*self.lam*n,q,.03,xtol=5e-18)
        rho=n*E-(E*E-q*q)**2/(2*self.lam)
        f=n/self.Ws
        S=self.Ws*(np.log1p(f)+f*np.log1p(1/f))
        return E,rho,S,f

    def pair_rhs(self,Y,closed):
        return super().rhs(np.r_[Y.ravel(),0.,0.],closed)

    def rhs_joint(self,y,closed,decay_scale):
        n=y[:-4]
        Y=n[:-1].reshape(-1,2)
        e,f,grad,_,_=self.weights(Y,closed)
        es,_,_,fs=self.scalar(n[-1])
        bounds=dict(coeff=e[self.start:],es=es)
        out=np.zeros_like(Y)
        scalar=0.;ent=0.
        for i,j,u,v,w in Tool.events(self.decay,bounds):
            i+=self.start;j+=self.start
            a=f[i,0]+f[i,1]*u
            b=f[j,0]+f[j,1]*v
            F=expit(a);B=expit(b)
            C=w*(F*(1-B)+fs*(F-B))*decay_scale
            n0=C.sum()
            out[i]-=np.r_[n0,C@u]
            out[j]+=np.r_[n0,C@v]
            scalar+=n0
            ent+=C@(a-b+np.log1p(1/fs))
        assert ent>=-1e-13
        work=(grad*out).sum()+es*scalar
        if closed: assert abs(work)<3e-13
        if self.pair_scale:
            z=self.pair_scale*self.pair_rhs(Y,closed)
            out+=z[:-2].reshape(-1,2);ent+=z[-2];work+=z[-1]
        assert abs(sum(out[:,0]))<3e-13*max(np.sum(abs(out[:,0])),1.)
        # Two independent ledgers separate total substitution and decay.
        return np.r_[out.ravel(),scalar,ent,work,scalar,out[self.r[0],0].sum()]

    def go_joint(self,closed,decay_scale=1.,duration=.14):
        init=np.r_[self.init.ravel(),self.Ns,0.,0.,0.,0.]
        def energy(y):
            E,S=self.weights(y[:-5].reshape(-1,2),closed)[-2:]
            _,r,s,_=self.scalar(y[-5])
            return E+r,S+s
        sol=solve_ivp(lambda _,y:self.rhs_joint(y,closed,decay_scale),(0.,duration),init,
                      method='DOP853',rtol=1e-10,atol=3e-13,max_step=.015)
        assert sol.success,sol.message
        end=sol.y[:,-1]
        e0,s0=energy(init);e1,s1=energy(end)
        assert abs(e1-e0-end[-3])<2e-12
        assert abs(s1-s0-end[-4])<3e-11
        if closed:assert abs(e1-e0)<2e-13
        dn=(end-init)[:-5].reshape(-1,2)[:,0]
        assert abs(dn.sum())<1e-13
        assert abs(end[-5]-self.Ns-end[-2])<1e-13
        assert abs(dn[self.r[0]].sum()-end[-1])<1e-13
        return dict(pairs=self.pair_scale,decay=decay_scale,closed=closed,duration=duration,
                    scalar_number_gain=end[-5]-self.Ns,energy_drift=e1-e0,integrated_work=end[-3],
                    entropy_gain=s1-s0,integrated_entropy=end[-4],combined_count_drift=dn.sum(),
                    species_number_change=[sum(dn[k]) for k in self.r])


def run(order):
    rows=[Joint(.005,order,p).go_joint(True,D) for p,D in ((0.,1.),(1.,0.),(.5,1.),(1.,1.),(2.,1.))]
    rows.append(Joint(.005,order,1.).go_joint(False,1.))
    return dict(scope='Common moment state, live Hartree scalar and reciprocal gauge energy; '
                      'incoming pair quadrature still just two test points and arbitrary clocks; no FRW.',
                order=order,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--order',type=int,default=10)
    a=p.parse_args()
    print(json.dumps(run(a.order),indent=2,allow_nan=False))
