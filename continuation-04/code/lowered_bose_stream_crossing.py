"""Independent monotone first crossing CDF for the narrow outgoing patch.
Valid before return; ignores exponentially tiny negative-mu and >6sigma tails.
"""
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from lowered_antonov_hessian import Halo

def cdf(model='linear',t=10,nr=128,nq=32):
 h=Halo(model=model);e0=h.Qt/h.rt;s=-.9*e0;rb=.8*h.rt;sig=.03*h.rt;mb=.8;sm=.05
 x,w=leggauss(nq);u=(x+1)/2;w=w/2
 x2,w2=leggauss(nr);r0=rb+6*sig*x2;mass=r0*r0*np.sqrt(h.psi(r0)-s)*np.exp(-.5*((r0-rb)/sig)**2)*w2
 vals=[]
 norm=erf((1-mb)/(np.sqrt(2)*sm))-erf(-mb/(np.sqrt(2)*sm))
 for r in r0:
  v0=2*(float(h.psi(r))-s);rr=r+(h.rt-r)*u*u;v=2*(h.psi(rr)-s);dr=2*(h.rt-r)*u*w
  def delay(mu):return float(np.sum(dr/np.sqrt(np.maximum(v-r*r*v0*(1-mu*mu)/rr**2,1e-300))))
  if delay(1)>t: vals.append(0.)
  elif delay(0)<=t: vals.append(1.)
  else:
   m=brentq(lambda mu:delay(mu)-t,0,1,xtol=1e-12);vals.append((erf((1-mb)/(np.sqrt(2)*sm))-erf((m-mb)/(np.sqrt(2)*sm)))/norm)
 return float(mass@np.asarray(vals)/mass.sum())
if __name__=='__main__':
 for n in [64,128,256,512]: print(n,cdf(nr=n))
