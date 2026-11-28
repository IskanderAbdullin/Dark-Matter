"""Sampled-angle reciprocal photon-pair evolution with live hard Bose moments."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_hard_bose_cell_reconstruction_reference import Cell as Gamma
from early_pole_live_pauli_gauge_moment_decay_reference import Cell as Fermion


class Pair:
    def __init__(self,g,n,edges=None):
        self.c=g*g/np.pi**4
        self.cells=[]
        self.r=[];state=[];self.shape=[]
        self.edge_sets=edges or ([.5,.92,1.12,1.55,2.0],
                                 [.45,.65,.85,1.05,1.25,1.7,2.2],
                                 [.45,.69,.8,1.03,1.3,1.7,2.2])
        self.background=np.zeros(3)
        for k,(edges,m) in enumerate(zip(self.edge_sets,(0.,.05,.003),strict=True)):
            begin=len(state)
            for l,h in zip(edges[:-1],edges[1:]):
                if k==0:
                    G=Gamma(l,h,n)
                    f=1/np.expm1(G.p+.3)
                    state.append(G.moments(f))
                    self.shape.append(((l+h)/2,(h-l)/2))
                else:
                    G=Fermion(l,h,m,n)
                    state.append(G.initial(-1.5-.4*k-G.ec,-G.de))
                    self.shape.append((G.ec,G.de))
                self.cells.append(G)
            self.r.append(slice(begin,len(state)))
        self.shape=np.asarray(self.shape)
        self.init=np.array(state)
        self.z,self.w=roots_legendre(n)
        self.incoming=[(.80,1.30,-.36,0,2),(.995,1.38,-.1,1,2)]

    def weights(self,Y,closed):
        F=[];D=[];J=[];S=0.
        for k,U in enumerate(Y):
            if k<self.r[0].stop:
                a,b,j,d,s=self.cells[k].fit(U)
                v=np.r_[a,b]
            else:
                v,j,d,s=self.cells[k].fit(U)
            F.append(v);D.append(d);J.append(j);S+=s
        J=np.array(J)
        j=self.background+np.array([sum(J[sl]) for sl in self.r])
        partners=np.r_[sum(j[1:]),j[0]+j[1:]]
        inc=np.empty(len(J))
        for k,v in enumerate(self.r):inc[v]=partners[k]
        H=(np.asarray(D)*inc[:,None])*self.c
        E=self.shape+H
        energy=(self.shape*Y).sum()+.5*self.c*(j@partners)
        return E if closed else self.shape,np.array(F),E,energy,S

    def rhs(self,z,closed):
        Y=z[:-2].reshape(-1,2)
        E,F,grad,_,_=self.weights(Y,closed)
        out=np.zeros_like(E);ent=0.
        for q1,q2,ct,p1,p2 in self.incoming:
            U=(np.array([q1,q2])-self.shape[[p1,p2],0])/self.shape[[p1,p2],1]
            EE=E[[p1,p2],0]+E[[p1,p2],1]*U
            X=F[[p1,p2],0]+F[[p1,p2],1]*U
            a1,a2=1/np.expm1(X)
            P=np.sqrt(q1*q1+q2*q2+2*q1*q2*ct)
            for span in self.r[1:]:
                for i in range(span.start,span.stop):
                    A=self.cells[i]
                    for j in range(span.start,span.stop):
                        B=self.cells[j]
                        es=EE.sum()
                        L=max(E[i,0]-E[i,1],es-E[j,0]-E[j,1])
                        H=min(E[i,0]+E[i,1],es-E[j,0]+E[j,1])
                        if H<=L:continue
                        inv=lambda k,x:self.cells[k].pfrombare(self.shape[k,0]+self.shape[k,1]*(x-E[k,0])/E[k,1])
                        p0,p1_=inv(i,L),inv(i,H)
                        ev=lambda p:E[i,0]+E[i,1]*A.basis(p)
                        pk=lambda p:inv(j,es-ev(p))
                        roots=[p0,p1_]
                        for f in (lambda p:p-pk(p)-P,lambda p:pk(p)-p-P,lambda p:P-p-pk(p)):
                            nodes=np.linspace(p0,p1_,9)
                            test=f(nodes)
                            for v in np.flatnonzero(test[:-1]*test[1:]<0):
                                roots.append(brentq(f,nodes[v],nodes[v+1],xtol=3e-14))
                        uu=np.sort(roots)
                        for l,h in zip(uu[:-1],uu[1:],strict=True):
                            mid=(l+h)/2;k=pk(mid)
                            if abs(mid-k)>P or P>mid+k:continue
                            p=mid+(h-l)*self.z/2;k=pk(p)
                            # Direction drawn from exact energy+triangle support.
                            u,v=A.basis(p),B.basis(k)
                            ea,eb=ev(p),es-ev(p)
                            va=E[j,1]/B.de*k/np.hypot(k,B.m)
                            w=(h-l)/2*self.w*p*k/(P*ea*eb*va)
                            f1,f2=F[i,0]+F[i,1]*u,F[j,0]+F[j,1]*v
                            P1,P2=expit(f1),expit(f2)
                            forward=a1*a2*(1-P1)*(1-P2)
                            bwd=(1+a1)*(1+a2)*P1*P2
                            r=w*(forward-bwd)
                            count=r.sum()
                            out[p1]-=count*np.r_[1.,U[0]]
                            out[p2]-=count*np.r_[1.,U[1]]
                            out[i]+=np.r_[count,r@u]
                            out[j]+=np.r_[count,r@v]
                            ent+=r@(-sum(X)-f1-f2)
        work=(grad*out).sum()
        assert ent>=-1e-14 and abs(out[:,0].sum())<2e-13
        if closed:assert abs(work)<2e-13
        return np.r_[out.ravel(),ent,work]

    def go(self,closed):
        a=self.init.ravel()
        sol=solve_ivp(lambda t,z:self.rhs(z,closed),(0.,.08),np.r_[a,0.,0.],
                      method='DOP853',rtol=1e-10,atol=3e-13,max_step=.012)
        assert sol.success
        n=sol.y[:-2,-1]
        _,_,_,e0,s0=self.weights(self.init,closed)
        _,_,_,e1,s1=self.weights(n.reshape(-1,2),closed)
        assert abs(e1-e0-sol.y[-1,-1])<2e-12
        assert abs(s1-s0-sol.y[-2,-1])<2e-11
        if closed:assert abs(e1-e0)<2e-13
        dn=(n-a).reshape(-1,2)[:,0]
        assert abs(dn.sum())<1e-13
        me=np.inf
        fmin,fmax=1.,0.
        for zz in sol.y[:-2].T:
            _,F,_,_,_=self.weights(zz.reshape(-1,2),closed)
            G=F[self.r[0]]
            me=min(me,np.min(G[:,0]-abs(G[:,1])))
            R=F[self.r[0].stop:]
            occ=expit(R[:,0,None]+R[:,1,None]*np.array([-1.,1.]))
            fmin=min(fmin,np.min(occ));fmax=max(fmax,np.max(occ))
        assert me>0. and 0.<fmin<=fmax<1.
        return dict(closed=closed,energy_drift=e1-e0,integrated_work=sol.y[-1,-1],
                    accumulated_entropy=s1-s0,integrated_entropy=sol.y[-2,-1],
                    min_bose_exponent=me,fermi_minmax=[fmin,fmax],
                    combined_count_drift=dn.sum(),
                    species_number_changes=[sum(dn[s]) for s in self.r])


def run():
    return dict(scope='Two fixed incoming hard angular samples; continuous outgoing two-body '
                      'geometry, no spin amplitude or physical clock.',
                rows=[dict(order=n,results=[Pair(.005,n).go(b) for b in (False,True)])
                      for n in (10,20)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
