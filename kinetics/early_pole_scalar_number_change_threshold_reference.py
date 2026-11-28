"""Four cold real scalars into two: vacuum threshold and frozen snapshots."""
import argparse
import gzip
import itertools
import json
from pathlib import Path

import numpy as np

TRIPLETS=np.array([(0,)+a for a in itertools.combinations(range(1,6),2)])
LAMBDA=7.65175e-8


def metric(p):
    return p[...,0]**2-np.sum(p[...,1:]**2,axis=-1)


def amplitude(signed,m=1.):
    three=signed[:,TRIPLETS,:].sum(axis=2)
    return np.sum(1/(metric(three)-m*m),axis=1)


def threshold():
    p=np.array([[[1.,0.,0.,0.]]*4+[[-2.,3.**.5,0.,0.],[-2.,-3.**.5,0.,0.]]])
    assert abs(amplitude(p)[0]+1.)<1e-14
    den=np.sort((metric(p[:,TRIPLETS,:].sum(axis=2))-1)[0])
    assert np.allclose(den,[-4.]*6+[8.]*4)
    event=np.sqrt(3)/(2*24*256*np.pi)
    return dict(amplitude_stripped_in_m1=-1.,internal_denominators_m2=den.tolist(),
                event_coefficient_times_lambda4_over_m8=event,
                net_number_loss_coefficient_times_lambda4_over_m8=2*event)


def directions(rng,shape):
    q=rng.normal(size=shape+(3,))
    return q/np.linalg.norm(q,axis=-1,keepdims=True)


def particles(rng,n,q,w):
    v=rng.choice(q,size=(n,4),p=w)
    p=directions(rng,(n,4))*v[:,:,None]
    E=np.hypot(v,1.)
    return np.concatenate((E[:,:,None],p),axis=2)


def exact_weight(spec,m,q,prob,rng,n,extreme):
    old=particles(rng,n,q,prob)
    total=old.sum(axis=1)
    cov=total[:,1:]/total[:,:1]
    s=metric(total)
    P0=np.sqrt(s)/2
    pCOM=directions(rng,(n,))*np.sqrt(s/4-1)[:,None]
    vc=np.sum(pCOM*cov,axis=1)
    cv=np.sum(cov*cov,axis=1)
    gamma=total[:,0]/np.sqrt(s)
    # (gamma-1)/v^2 = gamma^2/(gamma+1) remains finite.
    boost=gamma*gamma/(gamma+1)
    A=pCOM+(boost*vc+gamma*P0)[:,None]*cov
    B=-pCOM+(-boost*vc+gamma*P0)[:,None]*cov
    out0=np.c_[gamma*(P0+vc),A]
    out1=np.c_[gamma*(P0-vc),B]
    Q=np.stack((out0,out1),axis=1)
    p=np.concatenate((old,-Q),axis=1)
    assert np.max(abs(p.sum(axis=1)))<1e-13
    assert np.max(abs(metric(Q)-1))<2e-13
    theta0=np.linalg.norm(Q[:,:,1:],axis=2)*m
    assert theta0.min()>=extreme[0] and theta0.max()<=extreme[1]
    f=np.maximum(spec['endpoint_scalar_occupations'],0.)
    edges=np.asarray(spec['scalar_momentum_edges'])
    oi=np.searchsorted(edges,theta0,side='right')-1
    assert np.max(oi)<len(f)
    f_out=f[oi]
    ii=np.searchsorted(edges,np.linalg.norm(old[:,:,1:],axis=2)*m,side='right')-1
    f_in=f[ii]
    assert np.min(f_in)>0
    # Normalize against the all-at-rest cross section with no final factors.
    kinematic=amplitude(p)**2*np.sqrt((1.-4./s)/.75)/old[:,:,0].prod(axis=1)
    # Net F-R per product f_in. Separate terms are well-conditioned here
    # because the two outgoing modes are far from the stimulated IR.
    stim=np.prod(1+f_out,axis=1)
    inverse=np.prod(f_out,axis=1)*np.prod(1+1/f_in,axis=1)
    out=kinematic*(stim-inverse)
    return dict(samples=n,mean_kinematic_ratio=float(kinematic.mean()),
                standard_error_kinematic_ratio=float(kinematic.std()/np.sqrt(n)),
                mean_forward_factor=float((kinematic*stim).mean()),
                mean_reverse_factor=float((kinematic*inverse).mean()),
                mean_net_factor=float(out.mean()),standard_error_net_factor=float(out.std()/np.sqrt(n)),
                sampled_outgoing_q_interval=[float(theta0.min()),float(theta0.max())])


def one(s,M,seed,n):
    q=np.array(s['scalar_momentum_nodes'])
    edges=np.array(s['scalar_momentum_edges'])
    num=np.maximum(np.array(s['endpoint_scalar_signed_shell_number_over_T3']),0.)
    f=np.maximum(s['endpoint_scalar_occupations'],0.)
    r=s['checkpoints'][-1]; m=r['scalar_pole_over_T']
    cold=(edges[1:]<=1e-5)&(f>0)
    N=num.sum(); nc=num[cold].sum()
    cut=edges[1:][cold].max()
    # 4 future-timelike mass-m vectors give sqrt(s)>=4m. Use
    # E_total<=4 Eh and |P_total|<=4 q_cut for uniform outgoing bounds.
    de=2*cut*np.sqrt(1-m*m/(4*(m*m+cut*cut)))
    low=2*m-de
    high=2*np.hypot(m,cut)+de
    out=(np.sqrt((low-m)*(low+m)),np.sqrt((high-m)*(high+m)))
    assert out[0]>cut
    h=np.sqrt(np.pi**2*106.75/90)*(M/r['x'])/2.435e18
    coef=threshold()['net_number_loss_coefficient_times_lambda4_over_m8']
    nn=coef*LAMBDA**4*nc**4/(m**8*N*h)
    rr=exact_weight(s,m,q[cold]/m,num[cold]/nc,np.random.default_rng(seed),n,out)
    return dict(x=r['x'],mass_ratio=s['mass_ratio'],parent_mass_GeV=M,cold_number_fraction=nc/N,
                pole_over_T=m,analytic_outgoing_q_range=list(out),
                formal_vacuum_threshold_loss_over_reference_H=nn,
                decay_net_growth_per_lna=r['net_log_growth_rate'],
                formal_net_loss_with_final_factors_over_Href=nn*rr['mean_net_factor'],
                relativistic_phase_space_test=rr)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    p.add_argument('--samples',type=int,default=32768)
    p.add_argument('--parent-mass-GeV',type=float,default=86.)
    p.add_argument('--seed',type=int,default=50701)
    a=p.parse_args()
    out=dict(threshold=threshold())
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        out[f.name]=one(json.loads(raw),a.parent_mass_GeV,a.seed,a.samples)
    print(json.dumps(out,indent=2,allow_nan=False))
