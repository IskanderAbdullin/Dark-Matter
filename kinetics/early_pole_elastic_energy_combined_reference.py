"""Live-parent exact projected decay plus isolated parent mixing."""
import argparse
import json

import numpy as np

from early_pole_friedmann_energy_control_reference import FriedmannStencil
from early_pole_parent_elastic_relaxation_reference import ElasticBath


class ConservativeMix(FriedmannStencil):
    elastic_target=ElasticBath.elastic_target
    elastic_checks=ElasticBath.elastic_checks

    def __init__(self,ns,nf,nparent,rate,full_H):
        super().__init__(ns,nf,nparent,106.75,full_H)
        self.rate=rate

    def at(self,x,m,y,jac=False):
        v=super().at(x,m,y,jac)
        if self.rate:
            p=y[self.iP:]
            if jac:
                v[self.iP:,self.iP:]+=self.rate*(self.elastic_target(x,p,True)-np.eye(self.np))
            else:
                v[self.iP:]+=self.rate*(self.elastic_target(x,p)-p)
        return v

    def run(self,xend):
        checks=self.elastic_checks()
        out=super().run(xend)
        out['elastic_rate_over_radiation_H']=self.rate
        out['isolated_target_checks']=checks
        return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=65)
    p.add_argument('--fermion',type=int,default=90)
    p.add_argument('--parent',type=int,default=72)
    p.add_argument('--rate',type=float,default=100.)
    p.add_argument('--full-H',action='store_true')
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    assert a.rate>=0.
    d=ConservativeMix(a.scalar,a.fermion,a.parent,a.rate,a.full_H)
    print(json.dumps(d.run(a.xend),indent=2,allow_nan=False))
