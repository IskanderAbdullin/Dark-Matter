"""Closed finite-bath FD update with direct timelike and elastic heat quadratures."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from boost_axis_elastic_heat_reference import integrate as elastic, leptons
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_finite_bath_discrete_entropy import run as trajectories
from fd_portal_tabulated_bath_discrete_entropy import Tabulated
from sm_dark_entropy_exchange_reference import StandardModel


def elastic_heat(T,Td,eta,mesh,*,linear=False):
    n,cut=mesh
    return -sum(elastic(T,Td,eta,m,n=n,cut=cut,linear=linear) for m in leptons)


class Combined(Tabulated):
    def flux(self,p,u):
        y,L,g,T,A=self.point(u)
        a,s=super().flux(p,u)
        mesh=(18,29.) if self.dn==112 else (25,36.)
        q=elastic_heat(T,np.exp(u[1]),u[0],mesh)
        d=g[1]*q
        if d<0. or (d==0. and g[1]!=0.):
            raise ArithmeticError('Elastic entropy must have the direction of its thermal force')
        b=a+np.array([0.,q])
        if abs(g@b-(s+d))>8e-14*max(max(abs(b)),abs(s+d)):
            raise ArithmeticError('Timelike and elastic source entropy did not close')
        return b,s+d


def check_local(eos):
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,eta,r in ((.05,.3,1.16),(.3,-.75,1.16),(.6882487043,.9,.86)):
        u=np.array([eta,np.log(T*r)])
        a,b=[Combined(T,u,n,eos) for n in (112,168)]
        x,y=[c.flux(z,u) for c,z in ((a,p),(b,q))]
        base,entropy=Tabulated(T,u,168,eos).flux(q,u)
        if (y[0][0]!=base[0] or y[1]<entropy or np.max(abs(x[0]/y[0]-1.))>2e-5):
            raise ArithmeticError('Elastic change to number or full paired-quadrature mismatch')
        parts=[]
        for mesh in ((18,29.),(25,36.)):
            target=-elastic_heat(T,T,0.,mesh,linear=True)
            h=6e-5
            low=elastic_heat(T,T*(1-h),0.,mesh)
            high=elastic_heat(T,T*(1+h),0.,mesh)
            r0=(high-low)/(2*h)/target
            if abs(r0-1.)>3e-6:
                raise ArithmeticError('Small displacement does not recover elastic thermal derivative')
            parts.append(target*T)
        if abs(parts[0]/parts[1]-1)>2e-5:
            raise ArithmeticError('Two orders did not resolve positive elastic heat mobility')
        rows.append(dict(T_SM=T,initial_T_ratio=r,eta=eta,combined_number_and_heat=y[0].tolist(),
                         timelike_number_and_heat=base.tolist(),elastic_heat_mobility_at_common_state=parts[1],
                         elastic_entropy_fraction_at_initial_state=1.-entropy/y[1]))
    return rows


def run(path,local_only=False):
    eos=StandardModel(path)
    check=check_local(eos)
    trace=None if local_only else trajectories(lambda T,u,n:Combined(T,u,n,eos))['rows']
    return dict(phase_space_checks=check,combined_trajectories=trace,
                table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Local closed-system thermal FD benchmark with direct charged-lepton elastic and timelike vacuum terms, not a complete interaction history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('--local-only',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.local_only),indent=2,allow_nan=False))
