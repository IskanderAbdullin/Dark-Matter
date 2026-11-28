"""Late photon-temperature clock from the external central radiation table."""
import argparse,json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from sm_dark_entropy_exchange_reference import MPL


class TabulatedRadiation:
    def __init__(self,file,start):
        t,s,g=np.loadtxt(file).T
        keep=(t>0.)&(t<=1.1)
        t,s,g=[x[keep] for x in (t,s,g)]
        assert np.all(np.diff(t)>0.) and min(g)>0. and min(s)>0.
        u=np.log(t)
        self.sg=PchipInterpolator(u,np.log(s),extrapolate=False)
        self.rg=PchipInterpolator(u,np.log(g),extrapolate=False)
        self.lo,self.initial=np.log(t[0]),np.log(start)
        assert u[0]<self.initial<u[-1]
        self.s0=self.entropy_log(self.initial)

    def entropy_log(self,u):
        return float(3*u+self.sg(u))

    def state(self,b):
        assert b>=1.
        if b==1.: u=self.initial
        else:
            target=self.s0-3*np.log(b)
            if target<self.entropy_log(self.lo):
                raise ValueError("late clock outside positive-temperature EOS table")
            u=brentq(lambda x:self.entropy_log(x)-target,self.lo,self.initial,xtol=4e-14)
        T=np.exp(u)
        rho=np.pi**2/30*np.exp(4*u+self.rg(u))
        defect=(self.entropy_log(u)-self.s0)+3*np.log(b)
        assert abs(defect)<2e-12
        return float(T),float(rho)


def test(path,start=.01199971):
    bg=TabulatedRadiation(path,start)
    a=[]
    for b in (1.,3.,12.,24.,100.,346.,410.):
        T,r=bg.state(b)
        H=np.sqrt(r/3)/MPL
        lower=np.sqrt(np.pi**2/45)*(start/b)**2/MPL
        assert H>lower
        a.append(dict(b=b,photon_temperature_GeV=T,Hubble_radiation_GeV=H,
                      radiation_H_over_freely_redshifted_photon_only=H/lower))
    return a


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos')
    x=p.parse_args()
    print(json.dumps(dict(rows=test(x.eos)),indent=2,allow_nan=False))
