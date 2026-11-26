"""Initial quantum C22 creation above a local energy threshold, isotropic NR contact scattering."""
import math,json,argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from lowered_antonov_hessian import Halo

def rate(H,psi,excess=0,n=48):
 if psi<=excess or psi<=0:return (0.,0.)
 x,w=leggauss(n);x=(x+1)/2;w=w/2
 a=excess/psi+(1-excess/psi)*x[:,None,None];b0=1+excess/psi-a;b=b0+(1-b0)*x[None,:,None]
 q1=psi*a;q2=psi*b;S=q1+q2;P=psi+excess;A=S/2
 cm=np.sqrt(np.maximum(0,P*(S-P)/(q1*q2)));c=cm*(2*x[None,None,:]-1)
 B=np.sqrt(np.maximum(0,A*A-q1*q2*c*c));ds=np.maximum(0,(P*(S-P)-q1*q2*c*c)/(B+P-A));lo=psi+P-S;low=q1*q2*c*c/(A+B);hi=psi-low
 if H.model=='king':integ=np.exp(lo)*np.expm1(ds)
 else:
  integ=(1-H.c)*ds+np.log(-np.expm1(lo-H.eta))-np.log(-np.expm1(hi-H.eta))
  if H.model=='quadratic':integ-=H.d*ds*(lo+hi)/2
 qfac=(psi**3/4)*(1-excess/psi)*(1-b0)*np.sqrt(a*b)*H.F(psi-q1)*H.F(psi-q2)
 rel=np.sqrt(np.maximum(0,S-2*np.sqrt(q1*q2)*c));weights=w[:,None,None]*w[None,:,None]*w[None,None,:]*qfac*cm*rel
 return float(np.sum(weights*ds/B)),float(np.sum(weights*integ/B))

def integrated(H,n=48,nr=160):
 u,w=leggauss(nr);u=(u+1)/2;w/=2;r=H.rt*u*u;weights=w*2*H.rt*u*r*r;psi=H.psi(r);rows=[]
 for e in [0,H.Qt/H.rt]:
  vals=np.array([rate(H,p,e,n=n) for p in psi]);res=weights@vals
  mass=float(weights@H.I(psi)/H.Ic);assert abs(mass/H.Qt-1)<1e-5
  rows.append(dict(excess=e,Jclass=res[0],Jbose=res[1],enhancement=res[1]/res[0],mass_integral=mass))
 G=6.67430e-11;hbar=1.054571817e-34;evkg=1.78266192e-36;V=44100.;R=.6*3.085677581491367e19;rho=V*V/(4*math.pi*G*R*R);sigmam=.125;gyr=1e9*365.25*86400
 m=(rho*2*math.pi**2*hbar**3/(2**1.5*V**3*H.Ic))**.25
 pref=.5*rho*sigmam*math.sqrt(2)*V/(H.Ic**2*H.Qt)*gyr
 for d in rows:d['escape_created_per_particle_per_Gyr_bose']=d['Jbose']*pref;d['classical']=d['Jclass']*pref
 return dict(model=H.model,delta=H.delta,eta=H.eta,n=n,nr=nr,rt=H.rt,Qt=H.Qt,m_eV=m/evkg,rho0_kgm3=rho,results=rows)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--model',default='linear');ap.add_argument('--delta',type=float,default=0);ap.add_argument('--n',type=int,default=32);ap.add_argument('--nr',type=int,default=96);a=ap.parse_args();H=Halo(model=a.model,delta=a.delta);print(json.dumps(integrated(H,a.n,a.nr)))
