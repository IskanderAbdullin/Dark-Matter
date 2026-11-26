"""Evolving homogeneous isotropic C22 plus first-order l=2 angular perturbation.
Exact angular covariance on discrete energies; no nonlinear angular terms,
spatial streaming, gravitational potential update or condensate.
"""
import argparse,json
import numpy as np
from scipy.sparse.linalg import expm_multiply
from lowered_bose_quadrupole import occ
from lowered_bose_harmonic_c22 import matrices

class Network:
    def __init__(self,N,Emax,psi,model):
        self.N=N;self.de=de=Emax/N;self.E=E=(np.arange(N)+.5)*de;q=np.sqrt(E);self.D=q*de/2;self.M=self.D*E*E/5
        self.f0=occ(E,psi,model)[0];self.psi=psi
        allids=[];allw=[];allcov=[]
        for s in range(2*N-1):
            pair=np.arange(max(0,s-N+1),min(N-1,s//2)+1);other=s-pair
            aa,bb=np.triu_indices(len(pair),0);ids=np.array([pair[aa],other[aa],pair[bb],other[bb]]);ei=E[ids];qi=q[ids]
            w=de**3*np.minimum(qi[0],qi[2])/4*np.where(ids[0]==ids[1],1,2)*np.where(ids[2]==ids[3],1,2)/2*np.where(aa==bb,.5,1.)
            S=ei[0]+ei[1];delta=ei-ei[[1,0,3,2]]
            p0=np.maximum(qi[1]-qi[0],qi[3]-qi[2]);p1=np.minimum(qi[0]+qi[1],qi[2]+qi[3]);m2=(p0*p0+p0*p1+p1*p1)/3
            m4=(p0**4+p0**3*p1+p0*p0*p1*p1+p0*p1**3+p1**4)/5
            inv2=np.divide(1,p0*p1,out=np.zeros_like(p0),where=p0>0);inv4=np.divide(m2,p0**3*p1**3,out=np.zeros_like(p0),where=p0>0)
            ap2=(m2+2*delta+delta**2*inv2)/4;cov=np.empty((4,4,len(w)))
            for j in range(4):
                for i in range(4):
                    if i==j:dot2=ei[i]**2
                    elif i//2==j//2:dot2=(m4-2*S*m2+S*S)/4
                    else:
                        x,y=delta[i],delta[j];prod=(m4+2*(x+y)*m2+x*x+y*y+4*x*y+2*x*y*(x+y)*inv2+x*x*y*y*inv4)/16
                        dot2=1.5*prod+.5*(ei[i]*ei[j]-ei[i]*ap2[j]-ei[j]*ap2[i])
                    cov[j,i]=(3*dot2-ei[i]*ei[j])/10
            allids.append(ids.astype(np.int32));allw.append(w);allcov.append(cov)
        self.ids=np.concatenate(allids,axis=1);self.w=np.concatenate(allw);self.cov=np.concatenate(allcov,axis=2)
    def deriv(self,f,x):
        ids=self.ids; a,b,c,d=f[ids];w=self.w;N=self.N
        J=w*(a*b*(1+c+d)-c*d*(1+a+b));dJ=np.array([b*(1+c+d)-c*d,a*(1+c+d)-c*d,a*b-d*(1+a+b),a*b-c*(1+a+b)])
        y=np.einsum('jik,ik->jk',self.cov,dJ*x[ids],optimize=False)*w
        g=np.zeros(N);h=np.zeros(N)
        for j,sgn in enumerate([-1,-1,1,1]):
            g+=sgn*np.bincount(ids[j],weights=J,minlength=N);h+=sgn*np.bincount(ids[j],weights=y[j],minlength=N)
        return g/self.D,h/self.M
    def entropy(self,f):return self.D@((1+f)*np.log1p(f)-f*np.log(np.maximum(f,1e-300)))
    def report(self,t,f,x):
        b=self.E*x;lower=np.minimum(f+b,f-b/2);mask=f>1e-14
        return dict(t=t,N=float(self.D@f),E=float((self.D*self.E)@f),S=float(self.entropy(f)),Q=float(self.M@x),min_f=float(f.min()),min_angular_f=float(lower.min()),max_relative_angular=float(np.max(abs(b[mask])/f[mask])),outside=float(self.D[self.E>=self.psi]@f[self.E>=self.psi]),angular_outside=float(self.M[self.E>=self.psi]@x[self.E>=self.psi]))

def run(N,Emax,psi,model,tmax=1.,dtmax=.02,eps=.01):
    L=Network(N,Emax,psi,model);f=L.f0.copy();x=eps*f.copy();x0=x.copy();g,h=L.deriv(f,x);rate=-float(L.M@h/(L.M@x));Q0=float(L.M@x)
    E,D,_,_,Ks,_=matrices(N,Emax,psi,model,ls=(2,));mat=Ks[2]/L.M[:,None];check=float(np.max(abs(h-mat@x)))
    t=0.;rows=[L.report(t,f,x)];reject=0;steps=0
    for target in [.1,.5,1.,2.]:
        if target>tmax+1e-12:continue
        while t<target-1e-12:
            g,h=L.deriv(f,x);neg=g<0;lim=np.min(-f[neg]/g[neg]) if neg.any() else np.inf;dt=min(dtmax,target-t,.4*lim)
            while True:
                p=f+dt*g;z=x+dt*h;g2,h2=L.deriv(p,z);new=.5*(f+p+dt*g2);newx=.5*(x+z+dt*h2)
                if new.min()>=0 and L.entropy(new)>=L.entropy(f)-1e-12:break
                dt*=.5;reject+=1
            f,x=new,newx;t+=dt;steps+=1
        row=L.report(t,f,x);row.update(Q_frozen_rate=float(Q0*np.exp(-rate*t)),Q_frozen_operator=float(L.M@expm_multiply(t*mat,x0)));rows.append(row)
    return dict(N=N,Emax=Emax,psi=psi,model=model,eps=eps,dtmax=dtmax,initial_rate=rate,independent_matrix_action_max_abs=check,steps=steps,rejects=reject,rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=60);p.add_argument('--Emax',type=float,default=12);p.add_argument('--psi',type=float,default=3.9);p.add_argument('--model',default='linear');p.add_argument('--tmax',type=float,default=1);p.add_argument('--dtmax',type=float,default=.02);a=p.parse_args();print(json.dumps(run(a.N,a.Emax,a.psi,a.model,a.tmax,a.dtmax)))
