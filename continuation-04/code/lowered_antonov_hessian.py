"""Dimensionless energy-shell Antonov quadratic form for finite lowered Bose halos.
Newtonian isolated collisionless equilibria only. Not a spectral evolution solver.
"""
import math,json,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

class Halo:
 def __init__(self,eta=4,delta=0,model='linear',nq=192):
  self.eta=eta;self.psi0=eta-delta;self.delta=delta;self.model=model;self.c=1/np.expm1(eta);self.d=self.c*(1+self.c);self.u,self.w=leggauss(nq);self.u=(self.u+1)/2;self.w/=2
  if model=='king': self.psi0=eta-delta
  p=self.psi0;ps=np.unique(np.r_[np.linspace(0,p,1200),p-np.geomspace(1e-14,p,900),p]);vals=self.I(ps);self.Ic=float(vals[-1]);self.eos=PchipInterpolator(ps,vals/self.Ic)
  eps=1e-6
  def rhs(r,y):return [-y[1]/r**2,r*r*float(self.eos(np.clip(y[0],0,p)))]
  def edge(r,y):return y[0]
  edge.terminal=True;edge.direction=-1
  self.sol=solve_ivp(rhs,[eps,2000],[p-eps**2/6,eps**3/3],events=edge,dense_output=True,rtol=2e-10,atol=2e-12,max_step=.15,method='DOP853');assert self.sol.success and len(self.sol.t_events[0]);self.rt=float(self.sol.t_events[0][0]);self.Qt=float(self.sol.sol(self.rt)[1]);self.eps=eps
 def F(self,s):
  if self.model=='king':return np.expm1(s)
  v=self.eta-s; f=np.exp(v)*np.expm1(s)/(np.expm1(v)*np.expm1(self.eta))
  if self.model=='quadratic':
   small=self.d*(1+2*self.c)*s*s/2+self.d*(1+6*self.c+6*self.c**2)*s**3/6
   f=np.where(s<1e-5,small,f-self.d*s)
  return f
 def Fp(self,s):
  if self.model=='king':return np.exp(s)
  b=1/np.expm1(self.eta-s); v=b*(1+b)
  if self.model=='quadratic':v=np.where(s<1e-5,self.d*(1+2*self.c)*s+self.d*(1+6*self.c+6*self.c**2)*s*s/2,v-self.d)
  return v
 def integrals(self,p):
  p=np.asarray(p);q=np.sqrt(p[...,None])*self.u;y=self.eta-p;v=y[...,None]+q*q;s=p[...,None]-q*q
  if self.model=='king':f=self.F(s);return p**1.5*np.sum(self.w*self.u**2*f,axis=-1),np.sqrt(p)/2*np.sum(self.w*f,axis=-1)
  ysafe=np.maximum(y,1e-300);atan=np.arctan(np.sqrt(p/ysafe));reg=np.where(v<1e-3,-.5+v/12-v**3/720,1/np.expm1(v)-1/v)-self.c
  if self.model=='quadratic':reg-=self.d*s
  I=np.sqrt(p)-np.sqrt(y)*atan+p**1.5*np.sum(self.w*self.u**2*reg,axis=-1)
  W=.5*(atan/np.sqrt(ysafe)+np.sqrt(p)*np.sum(self.w*reg,axis=-1))
  f=self.F(s);I=np.where(p<.2,p**1.5*np.sum(self.w*self.u**2*f,axis=-1),I);W=np.where(p<.2,np.sqrt(p)/2*np.sum(self.w*f,axis=-1),W)
  return I,W
 def I(self,p):return self.integrals(p)[0]
 def W(self,r):return self.integrals(self.psi(r))[1]/self.Ic
 def psi(self,r):
  r=np.asarray(r);return np.maximum(0,np.where(r<self.eps,self.psi0-r*r/6,self.sol.sol(np.maximum(r,self.eps))[0]))
 def rs(self,s):
  if s<=0:return self.rt
  return brentq(lambda r:float(self.psi(r))-s,0,self.rt,xtol=1e-13)

def eval_fe(r,nodes):
 i=np.clip(np.searchsorted(nodes,r,side='right')-1,0,len(nodes)-2);t=(r-nodes[i])/(nodes[i+1]-nodes[i]);V=np.zeros((len(r),len(nodes)));j=np.arange(len(r));V[j,i]=1-t;V[j,i+1]=t;return V

def radial_nodes(rmax,n):return rmax*(np.linspace(0,1,n+1)**2)

def quad_segments(nodes,ng=8):
 u,w=leggauss(ng);dr=np.diff(nodes);r=nodes[:-1,None]+dr[:,None]*(u+1)/2;return r.ravel(),(dr[:,None]*w/2).ravel()

def forms(H,n=48,ne=128):
 nodes=radial_nodes(H.rt,n);r,w=quad_segments(nodes);V=eval_fe(r,nodes);dr=np.diff(nodes)
 A=np.zeros((n+1,n+1));coeff=(nodes[1:]**3-nodes[:-1]**3)/(3*dr*dr)
 for i,c in enumerate(coeff):A[i:i+2,i:i+2]+=c*np.array([[1,-1],[-1,1]])
 A[-1,-1]+=H.rt
 W=H.W(r);Blocal=V.T@(V*(w*r*r*W)[:,None]);M=V.T@(V*w[:,None]);
 u,ws=leggauss(ne);u=(u+1)/2;ws/=2;s=H.psi0*(1-u*u);ws=ws*2*H.psi0*u*H.Fp(s)/(2*H.Ic);B=np.zeros_like(A)
 for sk,wk in zip(s,ws):
  rs=H.rs(sk);rn=np.r_[nodes[nodes<rs],rs];rr,ww=quad_segments(rn,8);weight=ww*rr*rr*np.sqrt(np.maximum(H.psi(rr)-sk,0));X=eval_fe(rr,nodes);tot=weight.sum();X-=weight@X/tot;B+=wk*(X.T@(X*weight[:,None]))
 out={'eta':H.eta,'delta':H.delta,'model':H.model,'rt':H.rt,'Qt':H.Qt,'Ic':H.Ic,'n':n,'ne':ne,'radial_Bmax':float(eigh(B,A,eigvals_only=True,subset_by_index=[n,n])[0]),'max_r2W':float(np.max(r*r*W))}
 out['radial_naive_Bmax']=float(eigh(Blocal,A,eigvals_only=True,subset_by_index=[n,n])[0]);
 for l in [1,2,3]:
  Al=A+l*(l+1)*M;Al[-1,-1]+=l*H.rt
  v=eigh(Blocal[1:,1:],Al[1:,1:],eigvals_only=True,subset_by_index=[n-2,n-1]);out[f'l{l}_Bmax']=float(v[-1]);out[f'l{l}_Bnext']=float(v[-2])
 return out

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--model',default='linear');ap.add_argument('--delta',type=float,default=0);ap.add_argument('--eta',type=float,default=4);ap.add_argument('--n',type=int,default=48);ap.add_argument('--ne',type=int,default=128);a=ap.parse_args();h=Halo(a.eta,a.delta,a.model);print(json.dumps(forms(h,a.n,a.ne)))
