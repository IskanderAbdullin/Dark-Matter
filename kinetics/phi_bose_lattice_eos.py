"""Measure finite-time stress tensor, virial and spectral diagnostic on classical scalar checkpoints."""
import math,json,argparse
import numpy as np
from scipy.special import ellipj
from scipy.integrate import quad
from phi_bose_classical_lattice import Lattice

def moments(s):
 u=s.real();us=s.w*abs(s.u)**2;vs=s.w*abs(s.v)**2
 K=np.sum(vs)/2;G=np.sum(s.k2*us)/2;M=s.a*np.sum(us)/2;V=s.b*np.mean(u**4)/4
 uv=np.real(np.sum(s.w*s.u.conj()*s.v));E=K+G+M+V;P=K-G/3-M-V
 omega=np.sqrt(s.k2+s.a+3*s.b*np.sum(us));J=(vs+omega**2*us)/(2*omega);J[0,0,0]=0
 # Hartree proxy, not physical conserved number or free-particle PSD
 rms=math.sqrt(np.sum(s.k2*J)/np.sum(J));grad_uv=np.sum(s.k2*us*(s.k2>4))/(2*G)
 return np.array([K,G,M,V,E,P,uv,rms,grad_uv,float(np.max(abs(u)))])

def avg_homogeneous(s):
 p=s.p; period=s.period
 x2=quad(lambda z:ellipj(z,p)[1]**2,0,period,epsabs=1e-11)[0]/period
 E0=s.a/2+s.b/4
 return (1-s.a*x2/E0)/3

def run(s,periods=16,stride=8):
 assert s.spp%stride==0
 start=s.i*s.dt;v0=moments(s);v=v0.copy()/2;maxerr=abs(v0[4]/s.E0-1);last=v0;peak=v0[-1]
 N=periods*s.spp//stride
 for i in range(N):
  s.step(stride);last=moments(s);v+=last*(.5 if i==N-1 else 1);maxerr=max(maxerr,abs(last[4]/s.E0-1));peak=max(peak,last[-1])
 v/=N;K,G,M,V,E,P,uv,rms,grad,amp=v;elapsed=N*stride*s.dt;B=(last[6]-v0[6])/elapsed
 return dict(A=s.A,N=s.N,z0=start,z1=s.i*s.dt,periods=periods,w=P/E,w_hom=avg_homogeneous(s),
  w_boundary_corrected=(P-B/3)/E,trace_residual=(E-3*P-2*M+B)/E,virial_residual=(2*K-2*G-2*M-4*V-B)/E,
  K=K/E,G=G/E,M=M/E,V=V/E,E=E,uv_gradient_share=grad,hartree_k_rms=rms,peak_abs_u=peak,max_sampled_energy_error=maxerr)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--periods',type=int,default=16);ap.add_argument('--stride',type=int,default=8);ap.add_argument('--output');a=ap.parse_args();s=Lattice.load(a.input);print(json.dumps(run(s,a.periods,a.stride)));s.save(a.output) if a.output else None
