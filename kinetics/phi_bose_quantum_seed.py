"""Resonant-band instantaneous vacuum covariance, linear estimate and seeded classical runs."""
import math,json,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj,ellipk
from numpy.polynomial.legendre import leggauss
from phi_bose_classical_lattice import Lattice
from phi_bose_narrow_sweep import thermal,m

def monodromy(A, k2):
 p=A*A/(1+2*A*A);a=1-2*p;b=2*p;period=4*ellipk(p)
 def f(t,Y):
  u=ellipj(t,p)[1];f=k2+a+3*b*u*u
  a0,b0,c0,d0=Y
  return [c0,d0,-f*a0,-f*b0]
 return solve_ivp(f,[0,period],np.eye(2).ravel(),rtol=2e-10,atol=2e-12).y[:,-1].reshape(2,2)

def linear(A,lam,N=80,T=2e7):
 p=A*A/(1+2*A*A);a=1-2*p;b=2*p;period=4*ellipk(p);lo=3*(1+A*A)/(1+2*A*A);hi=lo+3*p*p/(2*math.sqrt(1-p+p*p)+2-p)
 x,w=leggauss(N);k=(math.sqrt(hi)+math.sqrt(lo))/2+x*(math.sqrt(hi)-math.sqrt(lo))/2;w=w*(math.sqrt(hi)-math.sqrt(lo))/2*k*k/(2*math.pi**2)
 C=lam/(6*b);om=np.sqrt(k*k+a+3*b);M=np.array([monodromy(A,q*q) for q in k]);V=np.zeros((N,2,2));V[:,0,0]=C/(2*om);V[:,1,1]=C*om/2
 E0=a/2+b/4;d0=np.sum(w*(V[:,1,1]+om**2*V[:,0,0]))/2
 initial=np.sqrt(np.sum(w*V[:,0,0]));result={'A':A,'lambda':lam,'T':T,'seed_rms':initial,'energy_vac_band_over_bg':d0/E0,'det_error':np.max(abs(np.linalg.det(M)-1)), 'H_over_Omega':thermal(T)[1]/(m*math.sqrt(1+2*A*A))}
 for n in range(1,1000):
  V=M@V@M.transpose(0,2,1);E=np.sum(w*(V[:,1,1]+om**2*V[:,0,0]))/2
  D=(E-d0)/E0
  for label,target in [('one_percent',.01),('ten_percent',.1)]:
   if D>target and label not in result:result[label]={'periods':n,'z':n*period,'Dlinear':D,'Hdt':n*period*result['H_over_Omega'],'sigma':np.sqrt(np.sum(w*V[:,0,0]))}
  if 'ten_percent' in result:break
 return result

def init(A,lam,N=48,L=32,spp=128,seed=42):
 s=Lattice(N,L,A,seed=seed,steps_per_period=spp);s.u[:]=0;s.v[:]=0;C=lam/(6*s.b);rng=np.random.default_rng(seed)
 for i in range(-N//4+1,N//4):
  for j in range(-N//4+1,N//4):
   for k in range(-N//4+1,N//4):
    if (i,j,k)<=(0,0,0):continue
    k2=(2*np.pi/L)**2*(i*i+j*j+k*k)
    if not s.gaplo<k2<s.gaphi:continue
    om=math.sqrt(k2+s.a+3*s.b)
    u=math.sqrt(C/(4*om*L**3))*(rng.normal()+1j*rng.normal());v=math.sqrt(C*om/(4*L**3))*(rng.normal()+1j*rng.normal())
    if k>=0:s.u[i%N,j%N,k]=u;s.v[i%N,j%N,k]=v
    if k<=0:s.u[-i%N,-j%N,-k]=u.conjugate();s.v[-i%N,-j%N,-k]=v.conjugate()
 s.u[0,0,0]=1;s.f=s.force();s.history=[];s.E0=s.measure()['energy'];s.record();return s

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--A',type=float,default=2.72085);ap.add_argument('--lam',type=float,default=4.698884e-13);ap.add_argument('--N',type=int,default=48);ap.add_argument('--L',type=float,default=32);ap.add_argument('--periods',type=int,default=0);ap.add_argument('--save');ap.add_argument('--linear',action='store_true');q=ap.parse_args()
 if q.linear:
  for N in [40,80]: print(json.dumps(linear(q.A,q.lam,N=N)))
 if q.periods:
  import os
  s=Lattice.load(q.save) if os.path.exists(q.save) else init(q.A,q.lam,N=q.N,L=q.L)
  for _ in range(q.periods):s.step(s.spp);print(json.dumps(s.record()),flush=True)
  s.save(q.save)
