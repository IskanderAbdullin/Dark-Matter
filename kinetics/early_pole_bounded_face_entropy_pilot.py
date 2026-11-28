"""Check the quantum entropy balance of geometric occupation events."""
import json
from math import fsum

import numpy as np
from scipy.special import xlogy

from early_pole_bounded_face_extents_pilot import step,ssp
from early_pole_double_charge_stencil_pilot import chart,event
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def ent(f,d):
    x,y=f[:d.ns],f[d.ns:]
    b=np.log1p(x)
    keep=x>0
    b[keep]+=x[keep]*np.log1p(1/x[keep])
    f=-xlogy(y,y)-xlogy(1-y,1-y)
    return np.r_[b,f]


def run():
    results=[]
    for grid,shape in SPECS.items():
        arc=load(grid)[0]
        x=arc['x_end']
        d=MassiveFermion(*shape,arc['mass_ratio'],arc['elastic_rate'],True)
        charge_mesh(d,1/4096,1/4096)
        u=np.concatenate([remap(a,arc['endpoint_'+s+'_occupations'],b)[0]
                          for a,b,s in zip(archives(arc),parts(d),LABELS,strict=True)])
        _,ra,rb,m=state_moments(d,x,u)
        H=np.sqrt(d.radiation/(d.radiation+ra+rb))
        ch,_=chart(d,x,m,4)
        s=u[:d.ns]
        y=u[d.ns:]
        e=np.r_[np.log1p(1/s),np.log1p(-y)-np.log(y)]
        val,fo,inv=event(ch,u,d)
        if np.min(fo)<=0 or np.min(inv)<=0:
            raise ArithmeticError('The entropy tangent needs interior occupation in every retained event')
        parts_rate=(fo-inv)*(np.log(fo)-np.log(inv))
        a,b=H*fsum(parts_rate),H*fsum(e*val)
        if np.min(parts_rate)<0 or a<=0 or abs(a/b-1.)>2e-13:
            raise ArithmeticError(('Event entropy does not close',grid,a,b))
        table=[]
        for alg in (step,ssp):
            for n in (16,64,256):
                f=u.copy()
                positive=[]
                seed=ent(f,d)
                cur=seed
                for i in range(n):
                    w,_=alg(ch,d,f,H*1e-6/n)
                    new=ent(w,d)
                    positive.append(fsum((new-cur)*d.w))
                    f=w
                    cur=new
                exact=fsum((cur-seed)*d.w)
                if min(positive)<0 or abs(exact-fsum(positive))>1e-14:
                    raise ArithmeticError('Sampled finite step did not preserve the entropy trend')
                table.append(dict(scheme=alg.__name__,steps=n,minimum_entropy_increment=min(positive),
                                  total_entropy_increase=exact))
        results.append(dict(grid=grid,entropy_production_event_and_ledger=[a,b],entropy_relative_defect=abs(a/b-1.),
                            scope='Entropy is of the held collision chart, not thermodynamic equivalence of two cosmologies.',
                            trials=table))
    return results


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
