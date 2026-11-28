"""Repeat the closed FD step with a thermodynamically consistent tabulated reservoir."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fd_portal_finite_bath_discrete_entropy import Bath, run
from sm_dark_entropy_exchange_reference import StandardModel


class Tabulated(Bath):
    def __init__(self,T,u,dn,eos):
        self.dn=dn
        self.eos=eos
        self.E=eos.state(T)[0]+self.dark(u)[0][1]
        self.lo,self.hi=np.exp([eos.lo+2e-14,eos.hi-2e-14])

    def point(self,u):
        y,s,A=self.dark(u)
        r=self.E-y[1]
        T=brentq(lambda v:self.eos.state(v)[0]-r,self.lo,self.hi,xtol=1e-14,rtol=1e-14)
        ev,_,sv,_=self.eos.state(T)
        if abs((y[1]+ev)/self.E-1.)>2e-12:
            raise ArithmeticError('Tabulated visible energy inversion failed')
        return y,-s-sv,np.array([u[0],1./T-np.exp(-u[1])]),T,A

    def visible_temperature_slope_denominator(self,T):
        return T*T*self.eos.state(T)[3]


def study(path):
    eso=StandardModel(path)
    z=run(lambda T,u,dn:Tabulated(T,u,dn,eso))
    return dict(rows=z['rows'],table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Fixed-energy local closed collision test using entropic tabulated visible EOS; not an expanding history or abundance calculation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(study(p.parse_args().eos),indent=2,allow_nan=False))
