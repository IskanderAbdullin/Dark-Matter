"""Find which destination mesh changes the local archived production clock."""
import json
from itertools import permutations,product
from math import fsum

import numpy as np

from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,X,archives,load,parts,remap,MassiveFermion


def effects(v):
    out=np.zeros(3)
    for idx in permutations(range(3)):
        st=['third']*3
        before=np.log(v[tuple(st)])
        for i in idx:
            st[i]='fine'
            after=np.log(v[tuple(st)])
            out[i]+=(after-before)/6
            before=after
    return out.tolist()


def run():
    targets=list(product(SPECS,repeat=3))
    result=[]
    for donor in SPECS:
        r=load(donor)[0]
        x=r['x_end']
        vals=[r['endpoint_'+k+'_occupations'] for k in LABELS]
        clocks={}
        out=[]
        for target in targets:
            shape=[SPECS[z][i] for i,z in enumerate(target)]
            d=MassiveFermion(*shape,r['mass_ratio'],r['elastic_rate'],True)
            remaps=[remap(b,u,a) for b,u,a in zip(archives(r),vals,parts(d),strict=True)]
            Y=np.concatenate([p[0] for p in remaps])
            mass=d.d.gap(x,Y[:d.ns])
            n=fsum(d.d.w*Y[:d.ns])
            sc,fb,pa=np.split(d.w*d.rhs(np.log(x/X),Y),(d.ns,d.iP))
            s=fsum(sc)
            if s<=0 or max(abs(s-fsum(fb)),abs(s+fsum(pa)))>3e-14:
                raise ArithmeticError('Production clock not resolved with number budget')
            clocks[target]=s/n
            out.append(dict(target_axes=target,net_log_clock=s/n,mass=mass,
                            flux_fraction_under_1e_5=fsum(sc[d.d.ed[1:]<=1e-5])/s))
        c=effects(clocks)
        if abs(sum(c)-np.log(clocks[('fine',)*3]/clocks[('third',)*3]))>3e-14:
            raise ArithmeticError('Grid axis audit failed to close')
        result.append(dict(source_archive=donor,axis_order=LABELS,log_effects=c,data=out))
    return result


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
