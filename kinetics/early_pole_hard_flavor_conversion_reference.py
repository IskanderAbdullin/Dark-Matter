"""Frozen hard number and energy flux of P/Pbar to b/bbar via a vacuum vector."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def integrate(level, bounds, finite_mass, reverse, fixtures, incident_edges=False):
    nr, ni, nc, nz=[(2,4,5,3),(3,6,7,4),(4,8,9,5),(1,3,4,3)][level]
    names=['fine','third']
    species=['parent','fermion']
    masses=[.05,.003] if finite_mass else [0.,0.]
    if reverse:
        species.reverse()
        masses.reverse()
    ma,mb=masses
    ma2,mb2=ma*ma,mb*mb
    profiles=[]
    for g in names:
        z=load(f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz')
        profiles.append([(np.asarray(z[f'endpoint_{k}_momentum_edges']),
                          np.r_[0.,z[f'endpoint_{k}_occupations'],0.]) for k in species])
    low,high=bounds
    if not 0.<low<high:
        raise ValueError("Choose a positive increasing radial window.")
    base=[.45,.55,.65,.8,1.,1.25,1.32,1.56,1.7,2.2,2.8,3.7,5.,6.5,8.,10.,14.]
    important=[x for a in profiles for x,w in zip(a[0][0],np.abs(np.diff(a[0][1])),strict=True)
               if 1.25<x<1.6 and w>.012]
    if incident_edges:
        important=[x for a in profiles for x in a[0][0] if low<x<high]
    knots=np.unique(np.r_[low,[x for x in base+important if low<x<high],high])
    p,dp=gauss_each(knots,nr)
    ep=np.hypot(p,ma)
    K=(1.-2*ma2)/2
    floor=(K*K-ep*ep*ma2)/(ep*np.sqrt(K*K-ma2*ma2)+p*K)
    z,w=roots_legendre(nr)
    l=np.maximum(knots[:-1][None,:],floor[:,None])
    half=np.maximum(0.,(knots[1:]-l)/2)
    b=l[...,None]+half[...,None]*(1+z)
    a=np.broadcast_to(p[:,None,None],b.shape)
    ea,eb=np.hypot(a,ma),np.hypot(b,ma)
    radial=(dp*p*p/np.pi**2)[:,None,None]*(b*b*half[...,None]*w/np.pi**2)
    ss=np.maximum(1.,2*(ma2+ea*eb-a*b))
    ds=np.maximum(0.,2*(ma2+ea*eb+a*b)-ss)
    z,w=roots_legendre(ni)
    S=ss[...,None]+ds[...,None]*(1+z)/2
    betain=np.sqrt(np.maximum(1-4*ma2/S,0.))
    rel=S*betain/(2*(ea*eb)[...,None])
    # Initial charge labels are fixed (2 spin states per leg). Since
    # forward integration begins with one chosen flavor, do not halve
    # the signed forward-minus-inverse rate again.
    W=(radial*ds/(8*a*b))[...,None]*w*rel
    keep=W>0
    s=S[keep]
    bi=betain[keep]
    bf=np.sqrt(1-4*mb2/s)
    a,b=[np.broadcast_to(x[...,None],W.shape)[keep] for x in (a,b)]
    W=W[keep]
    e1,e2=np.hypot(a,ma),np.hypot(b,ma)
    E=e1+e2
    boost=np.sqrt(np.maximum(E*E-s,0.))
    A=np.clip(np.divide(e1-e2,boost*bi,out=np.zeros_like(E),where=boost>0),-1.,1.)
    cz,cw=roots_legendre(nc)
    ph,pw=roots_legendre(nz)
    f=lambda q,k:k[1][np.searchsorted(k[0],q)]
    inside_edges=np.unique(np.r_[low,high,*[r[1][0][(r[1][0]>low)&(r[1][0]<high)] for r in profiles]])
    pe=np.hypot(inside_edges,mb)
    result=np.zeros((4 if fixtures else 2,5))
    norm=0.
    for c,w in zip(cz,cw,strict=True):
        weight=w*W*(bf/bi)*(1+(bi*bf*c)**2+4*(ma2+mb2)/s)/(32*np.pi*s)
        norm+=weight.sum()
        velocity=(E+boost*bf*A*c)/2
        span=np.maximum(boost*bf*np.sqrt((1-c*c)*(1-A*A))/2,1e-30)
        for k in range(0,len(E),224):
            h=slice(k,k+224)
            u=velocity[h]
            r=span[h]
            en=E[h]
            boundaries=np.c_[(pe-u[:,None])/r[:,None],(en[:,None]-pe-u[:,None])/r[:,None],
                              np.broadcast_to([-1.,1.],(len(en),2))]
            z=np.sort(np.arccos(np.clip(boundaries,-1.,1.)),axis=1)
            hp=np.diff(z,axis=1)/2
            i,j=np.nonzero(hp>1e-14)
            half=hp[i,j]
            mid=(z[i,j]+z[i,j+1])/2
            p3=u[i]+r[i]*np.cos(mid)
            p4=en[i]-p3
            v=(p3>=pe[0])&(p3<=pe[-1])&(p4>=pe[0])&(p4<=pe[-1])
            i,half,mid,p3,p4=[x[v] for x in (i,half,mid,p3,p4)]
            WW=2*half*weight[h][i]/np.pi
            for ix,(inc,out) in enumerate(profiles):
                A0=f(a[h][i],inc);B0=f(b[h][i],inc)
                C0=f(np.sqrt(np.maximum((p3-mb)*(p3+mb),0.)),out)
                D0=f(np.sqrt(np.maximum((p4-mb)*(p4+mb),0.)),out)
                fw=WW*A0*B0*(1-C0)*(1-D0)
                bw=WW*C0*D0*(1-A0)*(1-B0)
                flow=fw-bw
                result[ix]+=np.r_[flow.sum(),flow@en[i],np.sum(abs(flow)),fw.sum(),bw.sum()]
            if fixtures:
                E3=u[i,None]+r[i,None]*np.cos(mid[:,None]+half[:,None]*ph)
                E4=en[i,None]-E3
                A0=expit(.3-e1[h][i,None]);B0=expit(.3-e2[h][i,None])
                for q,chem in enumerate((.3,-.2),2):
                    C0=expit(chem-E3); D0=expit(chem-E4)
                    wt=WW[:,None]*pw/2
                    fw=wt*A0*B0*(1-C0)*(1-D0)
                    bw=wt*C0*D0*(1-A0)*(1-B0)
                    flow=fw-bw
                    result[q]+=np.r_[flow.sum(),np.sum(flow*en[i,None]),np.sum(abs(flow)),fw.sum(),bw.sum()]
    exact=np.sum(W*(bf/bi)*(2+2*(bi*bf)**2/3+8*(ma2+mb2)/s)/(32*np.pi*s))
    assert abs(norm/exact-1)<1e-13
    assert np.max(abs(result[:,0]-result[:,3]+result[:,4]))<2e-13
    if fixtures:
        assert max(abs(result[2,:2]))<1e-16 and min(result[3,:2])>0.
    H=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
    factor=.005**4/H
    return dict(scope='Strict finite vacuum s-channel with archived FD states; no soft/collinear or late rates.',
                order=[nr,ni,nc,nz],bounds=bounds,all_incident_edges=incident_edges,orientation=species,fermion_masses=masses,
                rows={key:dict(outgoing_species_count_gain=float(2*row[0]*factor),
                               outgoing_species_energy_gain=float(row[1]*factor),
                               absolute_event_flux=float(row[2]*factor),forward_event_flux=float(row[3]*factor),
                               reverse_event_flux=float(row[4]*factor))
                      for key,row in zip(names+(['common_T_mu','different_mu'] if fixtures else []),result,strict=True)})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,choices=range(4),default=0)
    p.add_argument('--bounds',nargs=2,type=float,default=[.45,8.])
    p.add_argument('--finite-mass',action='store_true')
    p.add_argument('--reverse',action='store_true')
    p.add_argument('--without-fixtures',dest='fixtures',action='store_false')
    p.add_argument("--incident-edges",action="store_true",help="split at every native incident occupation edge")
    a=p.parse_args()
    print(json.dumps(integrate(**vars(a)),indent=2,allow_nan=False))
