"""Fixed Minkowski classical real field: exact linear spectral rotation + dealiased quartic kick."""
import math,time,json,os,argparse
import numpy as np
from scipy.fft import rfftn,irfftn,fftfreq,rfftfreq
from scipy.special import ellipk,ellipj
from scipy.integrate import solve_ivp

class Lattice:
 def __init__(self,N=48,L=32,A=2.72085,seed=42,amp=1e-4,steps_per_period=64,workers=4):
  self.N=N;self.L=L;self.A=A;self.seed=seed;self.amp=amp;self.workers=workers
  self.p=A*A/(1+2*A*A);self.a=1/(1+2*A*A);self.b=2*A*A/(1+2*A*A)
  self.period=4*ellipk(self.p);self.dt=self.period/steps_per_period;self.spp=steps_per_period
  nx=fftfreq(N)*N;nz=rfftfreq(N)*N
  ii,jj,kk=np.meshgrid(nx,nx,nz,indexing='ij',sparse=True)
  self.k2=(2*np.pi/L)**2*(ii**2+jj**2+kk**2)
  # cubic de-aliasing and exact quartic quadrature: 4*K<N, retain strict half-Nyquist.
  self.mask=(abs(ii)<N/4)&(abs(jj)<N/4)&(abs(kk)<N/4)
  self.w=np.ones(N//2+1);self.w[1:-1]=2
  self.omega=np.sqrt(self.k2+self.a);self.cs=np.cos(self.omega*self.dt);self.sn=np.sin(self.omega*self.dt)
  self.gaplo=3*(1+A*A)/(1+2*A*A);self.gaphi=self.gaplo+3*self.p*self.p/(2*math.sqrt(1-self.p+self.p*self.p)+2-self.p)
  self.band=(self.k2>self.gaplo)&(self.k2<self.gaphi)
  rng=np.random.default_rng(seed);c=np.zeros((N,N,N),complex);kmax=math.ceil(math.sqrt(self.gaphi)*L/(2*np.pi))
  # Same integer Fourier coefficients at every N and Hermitian pairs.
  for i in range(-kmax,kmax+1):
   for j in range(-kmax,kmax+1):
    for k in range(-kmax,kmax+1):
     if (i,j,k)<=(0,0,0):continue
     k2=(2*np.pi/L)**2*(i*i+j*j+k*k)
     if self.gaplo<k2<self.gaphi:
      z=rng.normal()+1j*rng.normal();c[i%N,j%N,k%N]=z;c[-i%N,-j%N,-k%N]=z.conjugate()
  norm=np.sqrt(np.sum(abs(c)**2));assert norm>0
  c*=amp/norm;c[0,0,0]=1
  self.u=c[:,:,:N//2+1].copy();assert np.max(abs(self.u[~self.mask]))==0
  self.v=np.zeros_like(self.u);self.f=self.force();self.i=0;self.history=[];self.E0=self.measure()['energy'];self.record()
 def real(self):return irfftn(self.u,s=(self.N,)*3,norm='forward',workers=self.workers)
 def force(self):
  u=self.real();return self.b*self.mask*rfftn(u*u*u,norm='forward',workers=self.workers)
 def step(self,n):
  for _ in range(n):
   self.v-=.5*self.dt*self.f
   u=self.cs*self.u+(self.sn/self.omega)*self.v
   self.v=self.cs*self.v-(self.omega*self.sn)*self.u;self.u=u
   self.f=self.force();self.v-=.5*self.dt*self.f;self.i+=1
 def measure(self):
  u=self.real();us=abs(self.u)**2;vs=abs(self.v)**2
  E=.5*np.sum(self.w*(vs+(self.k2+self.a)*us))+self.b/4*np.mean(u**4)
  mean=float(self.u[0,0,0].real);vel=float(self.v[0,0,0].real)
  Emean=.5*vel*vel+.5*self.a*mean*mean+self.b/4*mean**4
  fluct=np.sum(self.w*us)-mean*mean;spec=self.w*(vs+self.k2*us)
  uv=float(np.sum(spec[self.k2>4])/max(np.sum(spec)-vel*vel,1e-300))
  return dict(step=self.i,z=self.i*self.dt,energy=float(E),mean=mean,sigma=math.sqrt(max(0,fluct)),nonmean=float(1-Emean/E),uv_k_gt_2=uv,band_variance=float(np.sum(self.w*us*self.band)))
 def record(self):
  d=self.measure();d['energy_rel']=d['energy']/self.E0-1;self.history.append(d);return d
 def save(self,path):
  np.savez(path,u=self.u,v=self.v,f=self.f,i=self.i,E0=self.E0,history=json.dumps(self.history),config=json.dumps(dict(N=self.N,L=self.L,A=self.A,seed=self.seed,amp=self.amp,steps_per_period=self.spp,workers=self.workers)))
 @classmethod
 def load(cls,path):
  with np.load(path) as z:
   s=cls(**json.loads(str(z['config'])));s.u=z['u'];s.v=z['v'];s.f=z['f'];s.i=int(z['i']);s.E0=float(z['E0']);s.history=json.loads(str(z['history']));return s

def main():
 p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=48);p.add_argument('--L',type=float,default=32);p.add_argument('--A',type=float,default=2.72085);p.add_argument('--spp',type=int,default=64);p.add_argument('--amp',type=float,default=1e-4);p.add_argument('--periods',type=int,default=4);p.add_argument('--save',required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 s=Lattice.load(a.save) if os.path.exists(a.save) else Lattice(a.N,a.L,a.A,amp=a.amp,steps_per_period=a.spp,workers=a.workers)
 t=time.monotonic()
 for _ in range(a.periods):s.step(s.spp);d=s.record();print(json.dumps(d),flush=True)
 s.save(a.save);print('elapsed',time.monotonic()-t,flush=True)
if __name__=='__main__':main()
