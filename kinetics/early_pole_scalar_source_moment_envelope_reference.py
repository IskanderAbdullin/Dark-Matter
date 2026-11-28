"""Exact frozen-pole source extrema at fixed scalar number and energy."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from scipy.optimize import linprog

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion


def envelope(q,K,S,mom):
    """Enumerate all at-most-two-point vertices at the specified mean K."""
    lo,hi=np.triu_indices(len(q),1)
    e0,e1=K[lo],K[hi]
    good=(e0<=mom)&(e1>=mom)
    lo,hi,e0,e1=(v[good] for v in (lo,hi,e0,e1))
    lam=(mom-e0)/(e1-e0)
    values=(1-lam)*S[lo]+lam*S[hi]
    return [(lo[j],hi[j],lam[j],values[j]) for j in (np.argmin(values),np.argmax(values))]


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
    for a,b,c,g in envelope(q[band],kinet[band],gain[band],mean):
        ix=np.flatnonzero(band)[[a,b]]
        u=np.where(band,0.,f)
        u[ix]=nt*np.array([1-c,c])/d.d.w[ix]
        z=np.r_[u,fd,P]
        expected=(nonband+nt*g)/total
        reevaluated=d.d.w@d.rhs(t,z)[:d.ns]/total
        pole=d.d.gap(x,u)
        assert abs((d.d.w@u)/total-1)<8e-15
        assert abs(pole/m-1)<3e-7
        extrema.append(dict(selected_q=q[ix].tolist(),shell_number_fractions=[1-c,c],
                            exact_frozen_rate_per_lna=float(expected),
                            recomputed_gap_rate_per_lna=float(reevaluated),
                            recomputed_gap_ratio=float(pole/m)))
    for j,sign in enumerate((1.,-1.)):
        parsed=linprog(sign*slope,A_eq=np.vstack((np.ones_like(keep),keep/mean)),b_eq=(1.,1.),
                       bounds=(0.,None),method='highs')
        assert parsed.success
        control=(nonband+nt*float(slope@parsed.x))/total
        assert abs(control-extrema[j]['exact_frozen_rate_per_lna'])<1e-8
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
