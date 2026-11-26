"""Signed initial quantum C22 flux across arbitrary energy thresholds via CM variables."""
import math,json,argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from lowered_antonov_hessian import Halo
from lowered_bose_collision_escape import rate

def occ(H,psi,q2):
 s=np.maximum(0,psi-q2);y=np.maximum(H.eta-psi+q2,1e-300)
 f=np.exp(y)*np.expm1(s)/(np.expm1(y)*np.expm1(H.eta))
 if H.model=='quadratic':f=np.where(s<1e-5,H.d*(1+2*H.c)*s*s/2+H.d*(1+6*H.c+6*H.c**2)*s**3/6,f-H.d*s)
 return f

def flux(H,psi,threshold,n=24,kind='number',logmap=False):
 # threshold is q_cut² (can exceed psi); positive = initial increase of below-threshold count.
 if psi<=0 or (kind=='number' and (threshold<=0 or threshold>=2*psi)):return 0.
 x,w=leggauss(n);x=(x+1)/2;w/=2
 # A=C²+d², b=2Cd. Separate A below/above threshold to smooth sign and turning points.
 edges=sorted(set([0.,psi]+[v for v in [threshold/2,threshold,(threshold+psi)/2] if 0<v<psi])) if kind=='number' else [0.,psi]
 total=0.
 gap=max(H.eta-psi,1e-15);off=max(gap+threshold,1e-12)
 for amin,amax in zip(edges[:-1],edges[1:]):
  if logmap:
   lA=np.log1p((amax-amin)/(amin+off));A=amin+(amin+off)*np.expm1(lA*x[:,None,None,None]);wa=w[:,None,None,None]*(A+off)*lA
  else:A=amin+(amax-amin)*x[:,None,None,None];wa=w[:,None,None,None]*(amax-amin)
  th0=.5*np.arcsin(np.minimum(1,np.abs(A-threshold)/A)) if kind=='number' else np.zeros_like(A)
  th=th0+(np.pi/4-th0)*x[None,:,None,None];wt=w[None,:,None,None]*(np.pi/4-th0)
  B=A*np.sin(2*th);uc=np.minimum(1,(psi-A)/B)
  t=np.abs(A-threshold)/B if kind=='number' else np.ones_like(B)
  uhi=np.minimum(uc,np.minimum(t,1))
  if logmap:
   du=np.maximum(gap+A-B*uhi,1e-15);lu=np.log1p(B*uhi/du);u=uhi-du*np.expm1(lu*x[None,None,:,None])/B;wu=w[None,None,:,None]*(du+B*(uhi-u))*lu/B
  else:u=uhi*x[None,None,:,None];wu=uhi*w[None,None,:,None]
  f1=occ(H,psi,A+B*u);f2=occ(H,psi,A-B*u);Pu=f1*f2;Su=f1+f2
  vlo=np.broadcast_to(np.minimum(t,1),u.shape) if kind=='number' else u
  for low,high in [(vlo,np.maximum(vlo,uc)),(np.maximum(vlo,uc),np.ones_like(vlo))]:
   high=np.minimum(high,1);width=np.maximum(0,high-low)
   if logmap:
    dv=np.maximum(gap+A-B*high,1e-15);lv=np.log1p(B*width/dv);v=high-dv*np.expm1(lv*x[None,None,None,:])/B;wv=w[None,None,None,:]*(dv+B*(high-v))*lv/B
   else:v=low+width*x[None,None,None,:];wv=width*w[None,None,None,:]
   f3=occ(H,psi,A+B*v);f4=occ(H,psi,A-B*v);Pv=f3*f4;Sv=f3+f4
   net=Pu*(1+Sv)-Pv*(1+Su)
   delta=np.sign(A-threshold) if kind=='number' else 2*B*B*(v*v-u*u)
   total+=np.sum(8*A**2.5*(np.cos(th)**2*np.sin(th)**3+np.sin(th)**2*np.cos(th)**3)*wa*wt*wu*wv*net*delta)
 return float(total)

def run(H,n=24,nr=48,logmap=False):
 x,w=leggauss(nr);x=(x+1)/2;w/=2;results=[]
 G=6.6743e-11;R=.6*3.085677581491367e19;V=44100.;rho=V*V/(4*math.pi*G*R*R);pref=.5*rho*.125*math.sqrt(2)*V/(H.Ic**2*H.Qt)*1e9*365.25*86400
 for d in [.2,.1,.03,.01,.003]:
  s0=H.eta-d
  if s0>=H.psi0:continue
  rt=H.rs(s0);r=rt*x;ww=w*rt*r*r;psi=H.psi(r)
  v=np.array([flux(H,p,p-s0,n,logmap=logmap) for p in psi]);J=ww@v
  q,wq=leggauss(96);q=(q+1)/2;wq/=2
  qmax=np.sqrt(psi-s0);count=ww@(qmax**3*np.sum(wq*q*q*occ(H,psi[:,None],qmax[:,None]**2*q*q),axis=1))/(H.Ic*H.Qt)
  results.append(dict(deep_gap=d,J=J,fraction=count,per_halo_Gyr=pref*J,per_deep_Gyr=pref*J/count))
 return results
def central_limit(H):
 # Pointwise limiting soft-mode C/f, including inverse reactions with unoccupied partner.
 assert H.delta==0 and H.model in ('linear','quadratic')
 eta=H.eta;c=H.c
 def f(x):return float(H.F(eta-x)) if 0<x<eta else 0.
 def inn(S):
  if H.model=='linear' and S<=eta:return c*(1+2/np.expm1(S))-c*c
  z0=max(0.,1-eta/S);z1=min(1.,eta/S);f2=f(S)
  def fun(z):
   f3=f(S*z);f4=f(S*(1-z));return f3*f4-f2*(1+f3+f4)
  return quad(fun,z0,z1,epsabs=1e-10,epsrel=1e-8)[0]
 j1=quad(lambda S:S*.5*inn(S),0,eta,epsabs=1e-9)[0]
 j2=quad(lambda S:S*.5*inn(S),eta,2*eta,epsabs=1e-9)[0]
 rho=44100**2/(4*math.pi*6.6743e-11*(.6*3.085677581491367e19)**2)
 pref=rho*.125*math.sqrt(2)*44100/H.Ic*1e9*365.25*86400
 return dict(j_bound=j1,j_above_cutoff=j2,soft_fractional_rate_Gyr=pref*(j1+j2),deep_fraction_over_gap_squared=6**1.5*math.pi/(64*H.Ic*H.Qt))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--model',default='linear',choices=['linear','quadratic']);ap.add_argument('--delta',type=float,default=0);ap.add_argument('--n',type=int,default=24);ap.add_argument('--nr',type=int,default=40);ap.add_argument('--log-map',action='store_true');a=ap.parse_args();H=Halo(model=a.model,delta=a.delta)
 for p in [1,3.9]:
  f=flux(H,p,p,a.n,logmap=a.log_map);print('escape',p,-f,rate(H,p,0,48)[1],flush=True)
 print(json.dumps(run(H,a.n,a.nr,logmap=a.log_map)),flush=True)
 if a.delta==0 and a.model in ('linear','quadratic'):print(json.dumps(central_limit(H)),flush=True)
