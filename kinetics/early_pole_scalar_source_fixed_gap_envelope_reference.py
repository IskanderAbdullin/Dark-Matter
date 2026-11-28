"""Source extrema preserving scalar number, energy and inverse-energy gap."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from scipy.optimize import linprog

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion


def one(spec,lo,hi):
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=MassiveFermion(*shape[tuple(spec['cells'])],spec['mass_ratio'],0.,True)
    x=spec['x_end']
    t=np.log(x/X)
    f=np.maximum(np.asarray(spec['endpoint_scalar_occupations']),0.)
    fd=np.asarray(spec['endpoint_fermion_occupations'])
    P=np.asarray(spec['endpoint_parent_occupations'])
    q=d.d.p
    y=np.r_[f,fd,P]
    _E,D,pp,m=state_moments(d,x,y)
    e,ep,i,j,k,l,a,b,W,miss=d.channels(x,m)
    j=j-d.ns
    F,B=d.parent_products(y[k],y[l],a,b,False)
    do=F*(1-fd[j])-B*fd[j]
    scale=np.sqrt(d.radiation/(d.radiation+D+pp))
    seed=scale*np.sum(W*F*(1-fd[j]))
    gain=scale*np.bincount(i,weights=W*do,minlength=d.ns)/d.d.w
    num=d.d.w*f
    total=num.sum()
    kinet=q*q/(e+m)
    band=(d.d.ed[:-1]>=lo)&(d.d.ed[1:]<=hi)
    nt=num[band].sum()
    mean=kinet[band]@num[band]/nt
    nonband=seed+gain[~band]@num[~band]
    extrema=[]
    keep=kinet[band]
    slope=gain[band]
    height=keep**2/e[band]
    hh=height@num[band]/nt
    for sign in (1.,-1.):
        parsed=linprog(sign*slope,A_eq=np.vstack((np.ones_like(keep),keep/mean,height/hh)),
                       b_eq=(1.,1.,1.),bounds=(0.,None),method='highs')
        assert parsed.success,parsed.message
        mask=parsed.x>0
        r=parsed.x[mask]
        ix=np.flatnonzero(band)[mask]
        u=np.where(band,0.,f)
        u[ix]=nt*r/d.d.w[ix]
        z=np.r_[u,fd,P]
        expected=(nonband+nt*float(slope@parsed.x))/total
        reevaluated=d.d.w@d.rhs(t,z)[:d.ns]/total
        pole=d.d.gap(x,u)
        assert abs((d.d.w@u)/total-1)<2e-10
        assert abs(pole/m-1)<4e-10
        assert abs(expected-reevaluated)<2e-6*max(1,abs(expected))
        extrema.append(dict(selected_q=q[ix].tolist(),shell_number_fractions=r.tolist(),
                            exact_frozen_rate_per_lna=float(expected),
                            recomputed_gap_rate_per_lna=float(reevaluated),
                            recomputed_gap_ratio=float(pole/m)))
    old=(seed+gain@num)/total
    direct=d.d.w@d.rhs(t,y)[:d.ns]/total
    assert abs(old-direct)<1e-10*max(1.,abs(direct))
    assert extrema[0]['exact_frozen_rate_per_lna']<=old<=extrema[-1]['exact_frozen_rate_per_lna']
    return dict(x=x,mass_ratio=spec['mass_ratio'],cells=spec['cells'],
                redistributed_band_shell_count=int(band.sum()),redistributed_number_fraction=nt/total,
                initial_net_source_per_lna=float(old), frozen_pole=m,
                extreme_min=extrema[0],extreme_max=extrema[1])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    p.add_argument('--lower',type=float,default=1e-6)
    p.add_argument('--upper',type=float,default=1e-5)
    a=p.parse_args()
    out={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        out[f.name]=one(json.loads(raw),a.lower,a.upper)
    print(json.dumps(out,indent=2,allow_nan=False))
