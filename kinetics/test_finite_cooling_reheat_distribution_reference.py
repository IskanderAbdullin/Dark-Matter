"""Independent dense clipping tests for sparse decay and streaming bounds."""
import argparse
import json

import numpy as np

from finite_cooling_reheat_distribution_reference import (
    Background, SourceFold, TF, check_partial_angles, decay_bins, positive_range_sum,
    stream_bounds,
)
from three_species_decay_onsager_reference import M, P, Es, m, mp


def dense(edge,mass,c,w,g,red):
    cut=np.hypot(mass,edge/red)
    cos=np.clip((cut[None,:]-c[:,None])/w[:,None],-1.,1.)
    x=np.diff(cos,axis=1)
    lo=c[:,None]+w[:,None]*cos[:,:-1]
    width=x*w[:,None]
    weight=g[:,None]*x/2
    p2=np.maximum(0.,(lo-mass)*(lo+mass)+lo*width+width*width/3)*red**2
    return np.array([weight.sum(axis=0),(weight*p2).sum(axis=0)])


def check(table):
    K=SourceFold(table,(8,23,17),121,4)
    edge=np.r_[0.,np.geomspace(1e-5,110,81),np.inf]
    err=[]
    for T,red in ((.001,1.),(.003,.3),(.005,.15),(.008,.1),(.01,.07),(.025,.01),(.032,.006)):
        E,W=K.C[1].particles(T,28)
        width=np.sqrt((E-M)*(E+M))*P/M
        for mass,E0 in ((mp,Es),(m,M-Es)):
            c=E*(E0/M)
            v=decay_bins(edge,mass,c,width,W,red/TF)
            v0=dense(edge,mass,c,width,W,red/TF)
            zz=np.sum(abs(v-v0),axis=-1)/v0.sum(axis=-1)
            err.append(zz)
    assert np.max(err)<1e-10

    rng=np.random.default_rng(428)
    ni=258
    bounds=np.sort(rng.integers(0,78,(2,ni)),axis=0)
    z=np.exp(rng.uniform(-45,1,ni))
    fast=positive_range_sum(*bounds,z,77)
    values=np.arange(77)[:,None]
    ref=((values>=bounds[0])&(values<bounds[1]))@z
    rel=np.max(abs(fast/ref-1.))
    assert rel<1e-14
    # This synthetic distribution puts two samples in each of several bins.
    points=np.array([.025,.042,.08,.144,.36,.72,1.28,2.1,5.25,7.88])
    cut=np.array([0.,.07,.17,1.,3.,10.,np.inf])
    mass=mp
    bag=np.arange(1,len(points)+1,dtype=float)
    n=np.histogram(points,bins=cut,weights=bag)[0]
    p2=np.histogram(points,bins=cut,weights=bag*points**2)[0]
    bg=Background(table)
    bounds=stream_bounds(bg,cut,n,p2,mass,8)
    ai=bg.scale(TF)
    v=bg.streams(points*TF*ai,mass,ai,bg.aeq,8)
    mu=bag@v/bag.sum()
    assert bounds["mean_Mpc"][0]<mu<bounds["mean_Mpc"][1]
    index=np.searchsorted(np.cumsum(bag)/bag.sum(),[.1,.5,.9])
    lo,hi=np.array(bounds["quantile_brackets_Mpc"]).T
    assert ((lo<=v[index])&(v[index]<=hi)).all()
    return dict(max_source_node_count_and_p2_l1_relative_difference=np.max(err,axis=0).tolist(),
                max_positive_range_sum_relative_error=float(rel),
                analytic_partial_angle_absolute_difference=check_partial_angles(),
                synthetic_exact_stream_mean_Mpc=float(mu),synthetic_stream_bounds=bounds)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    a=p.parse_args()
    print(json.dumps(check(a.eos),indent=2,allow_nan=False))
