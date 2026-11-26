"""Axisymmetric UU using continuous collision geometry and conservative
joint rounding of axial momenta. Mass-lumped hat quadrature."""
import argparse,json,numpy as np
from scipy.stats import qmc
from numpy.polynomial.legendre import leggauss
from lowered_bose_quadrupole import occ
from lowered_bose_axis_lattice import Network as BaseNetwork
from lowered_bose_axis_balance import balance

def hats(q,dp,K):
    x,w=leggauss(4);rows=[]
    for n in range(-K,K+1):
        p=n*dp;edges=np.unique(np.clip([p-dp,p,p+dp],-q,q));r=np.zeros(3)
        for lo,hi in zip(edges[:-1],edges[1:]):
            z=(hi+lo)/2+(hi-lo)*x/2;wt=w*(hi-lo)/2*np.maximum(1-abs(z-p)/dp,0);r+=np.array([wt.sum(),wt@z,wt@(z*z)])
        rows.append(r)
    return np.array(rows)

def pairprob(x,y,T):
    a=np.floor(x).astype(int);b=np.floor(y).astype(int);f=x-a;g=y-b;k=T-a-b;u=f+g
    p=np.where(u<=1,np.divide(f,u,out=np.zeros_like(f),where=u>0),np.divide(1-g,2-u,out=np.ones_like(f),where=u<2))
    p=np.where(k==0,0.,np.where(k==2,1.,p));return a,np.clip(p,0,1),k

def rounded_distribution(E,dp,K,samples):
    q=np.sqrt(E);a=max(abs(q[0]-q[1]),abs(q[2]-q[3]));b=min(q[0]+q[1],q[2]+q[3]);P=a+(b-a)*samples[:,0];mu=2*samples[:,1]-1;z=np.sqrt(1-mu*mu);S=P*mu/dp
    A=(P*P+E[0]-E[1])/(2*P);B=np.sqrt(np.maximum(0,E[0]-A*A));x=(A*mu+B*z*np.cos(2*np.pi*samples[:,2]))/dp;y=S-x
    A=(P*P+E[2]-E[3])/(2*P);B=np.sqrt(np.maximum(0,E[2]-A*A));v=(A*mu+B*z*np.cos(2*np.pi*samples[:,3]))/dp;w=S-v
    t0=np.floor(S).astype(int);frac=S-t0;M=2*K+1;ans=np.zeros(M**3);N=len(S)
    for T,wt in [(t0,1-frac),(t0+1,frac)]:
        n,p,k=pairprob(x,y,T);m,r,l=pairprob(v,w,T)
        assert np.all(((k>=0)&(k<=2))|(wt<1e-12)) and np.all(((l>=0)&(l<=2))|(wt<1e-12))
        for j in [0,1]:
            for i in [0,1]:
                n1=n+j;n2=T-n1;n3=m+i;n4=T-n3;weight=wt*(p if j else 1-p)*(r if i else 1-r)/N
                ok=weight>0;n1,n2,n3,n4,weight=[a[ok] for a in (n1,n2,n3,n4,weight)];assert np.all(abs(np.array([n1,n2,n3,n4]))<=K)
                key=((n1+K)*M+n2+K)*M+n3+K;rev=((K-n1)*M+K-n2)*M+K-n3
                ans+=np.bincount(key,weights=weight/2,minlength=M**3)+np.bincount(rev,weights=weight/2,minlength=M**3)
    key=np.flatnonzero(ans);n3=key%M-K;key2=key//M;n2=key2%M-K;n1=key2//M-K
    return np.array([n1,n2,n3,n1+n2-n3]),ans[key]

class Network(BaseNetwork):
    def __init__(self,NE,K,Emax,m=8,seed=0,bal=True,maxit=40):
        self.NE=NE;self.K=K;self.de=de=Emax/NE;self.dp=dp=np.sqrt(Emax)/K;self.Emax=Emax;en=(np.arange(NE)+.5)*de
        table=np.array([hats(np.sqrt(e),dp,K) for e in en]);ie,ip=np.where(table[:,:,0]>1e-16);self.ie=ie;self.ip=ip-K;self.N=N=len(ie);self.E=en[ie];self.p=self.ip*dp;mom=table[ie,ip];self.D=de/4*mom[:,0];self.Ebar=self.E;self.p2bar=mom[:,2]/mom[:,0];self.H=1.5*self.p2bar-.5*self.E;self.physical_p=mom[:,1]/mom[:,0]
        tab=np.full((NE,2*K+1),-1,dtype=int);tab[ie,ip]=np.arange(N);el=[np.flatnonzero(ie==e) for e in range(NE)];targets=np.zeros((NE,N))
        for e,ids in enumerate(el):targets[e,ids]=self.D[ids]/self.D[ids].sum()
        samples=qmc.Sobol(4,scramble=True,seed=seed).random_base2(m);self.groups=0;self.missing=0;self.nonconverged=0;self.maxerr=0.;self.meanerr=0.;self.meanL1change=0.;allids=[];allw=[]
        for s in range(2*NE-1):
            pairs=np.arange(max(0,s-NE+1),min(NE-1,s//2)+1);other=s-pairs
            for a,b in zip(*np.triu_indices(len(pairs),0)):
                es=np.array([pairs[a],other[a],pairs[b],other[b]]);E=en[es];momidx,v=rounded_distribution(E,dp,K,samples);ids=tab[es[:,None],momidx+K];assert np.all(ids>=0);raw=v/v.sum();v=raw
                if bal:
                    v,err,it,possible=balance(v,ids,targets[es],N,maxit=maxit);self.missing+=not possible;self.nonconverged+=err>1e-6;self.maxerr=max(self.maxerr,err);self.meanerr+=err;self.meanL1change+=float(np.sum(abs(v-raw)))
                self.groups+=1;W=de**3*np.sqrt(E).min()/8*(1 if es[0]==es[1] else 2)*(1 if es[2]==es[3] else 2)*(.5 if a==b else 1.)
                keep=v>0;allids.append(ids[:,keep].astype(np.int32));allw.append(W*v[keep])
        self.ids=np.concatenate(allids,axis=1);self.w=np.concatenate(allw);self.meanerr/=self.groups;self.meanL1change/=self.groups
    def isotropy(self,g):
        bar=np.bincount(self.ie,weights=self.D*g,minlength=self.NE)/np.bincount(self.ie,weights=self.D,minlength=self.NE)
        return float(self.D@abs(g-bar[self.ie])/max(self.D@abs(g),1e-300))

def run(NE,K,Emax,psi,model,m,seed,bal,maxit,tmax=.1):
    L=Network(NE,K,Emax,m,seed,bal,maxit);E=L.E;p=L.p;eps=.01;f0=lambda e:np.exp(-e) if model=='MB' else occ(e,psi,model)[0];f=f0(E)
    a=f0(E+np.expm1(-eps)*L.Ebar+(np.exp(2*eps)-np.exp(-eps))*L.p2bar);b=f0(E+np.expm1(eps)*L.Ebar+(np.exp(-2*eps)-np.exp(eps))*L.p2bar);bose=model!='MB';g=L.deriv(f,bose);ga=L.deriv(a,bose);gb=L.deriv(b,bose);rate=-float((L.D*L.H)@(ga-gb)/((L.D*L.H)@(a-b)))
    be=1/np.expm1(.5+E-.1*p);out=dict(NE=NE,K=K,Emax=Emax,model=model,bal=bal,m=m,seed=seed,states=L.N,events=len(L.w),groups=L.groups,missing_groups=L.missing,nonconverged_groups=L.nonconverged,max_marginal_L1=L.maxerr,mean_marginal_L1=L.meanerr,mean_weight_L1_change=L.meanL1change,rate=rate,isotropy=L.isotropy(g),Qdot=float((L.D*L.H)@g),outside_derivative=float(L.D[E>=psi]@g[E>=psi]),eq_res=float(abs(L.deriv(be)).max()),ledger=[float(L.D@g),float((L.D*E)@g),float((L.D*p)@g)])
    if not bose:return out
    f=a;out['initial']=L.report(f);t=0;steps=0;rejects=0
    while t<tmax-1e-12:
        g=L.deriv(f);neg=g<0;lim=np.min(-f[neg]/g[neg]) if neg.any() else np.inf;d=min(.02,tmax-t,.4*lim)
        while True:
            x=f+d*g;y=.5*(f+x+d*L.deriv(x))
            if x.min()>=0 and y.min()>=0 and L.report(y)['S']>=L.report(f)['S']-1e-12:break
            d*=.5;rejects+=1
        f=y;t+=d;steps+=1
    out.update(final=L.report(f),steps=steps,rejects=rejects,t=t,outside_final=float(L.D[E>=psi]@f[E>=psi]));return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=16);p.add_argument('--K',type=int,default=8);p.add_argument('--Emax',type=float,default=12);p.add_argument('--psi',type=float,default=3.9);p.add_argument('--model',default='linear');p.add_argument('--m',type=int,default=8);p.add_argument('--seed',type=int,default=0);p.add_argument('--maxit',type=int,default=40);p.add_argument('--raw',action='store_true');p.add_argument('--tmax',type=float,default=.1);a=p.parse_args();print(json.dumps(run(a.NE,a.K,a.Emax,a.psi,a.model,a.m,a.seed,not a.raw,a.maxit,a.tmax)))
