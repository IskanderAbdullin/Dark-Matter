"""Gross finite-width visible production and stable reverse loss on archives."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import cropped,snapshot_functions
from empty_lepton_pair_spectrum_reference import Source
from late_radiation_background_reference import TabulatedRadiation, MPL
from three_species_decay_onsager_reference import m
from timelike_portal_onsager_coverage_reference import V,width


def rates(T,s,order,npolar,tail):
    edges=np.array(s["lab_momentum_edges_GeV"])
    _,_,chi,eg=snapshot_functions(s)
    source=Source(order,dark_mass=m,tail=tail)
    q,p,v,N,A=source.event_moments(T,include_boltzmann=False)
    # Geometry in each boost block covers the entire represented reverse
    # input energy range when the excess cutoff reaches its upper edge.
    assert 2*np.hypot(m,edges[-1])<2*m+T*tail
    M2=q*q-p*p
    z,w=roots_legendre(npolar)
    accum=np.zeros((4,2))
    q=q.ravel(); N=N.ravel();A=A.ravel()
    h=(p*np.sqrt(v)[:,None]/2).ravel()
    M=M2.ravel()**.5
    for i in range(0,len(q),192):
        sl=slice(i,i+192)
        tot=q[sl,None]
        center=tot/2.
        dh=h[sl,None]
        bands=np.sort(np.clip(np.concatenate([-np.ones_like(center),(eg-center)/dh,
                                 (center-eg)/dh,np.ones_like(center)],axis=1),-1.,1.),axis=1)
        r=(bands[:,1:]+bands[:,:-1])/2
        d=(bands[:,1:]-bands[:,:-1])/2
        x=(r[:,:,None]+d[:,:,None]*z).reshape(len(r),-1)
        weight=(d[:,:,None]*w).reshape(len(r),-1)
        W=weight*(N[sl,None]+A[sl,None]*(3*x*x-1)/2)
        assert np.min(W)>=0
        a=chi(center+dh*x); b=chi(center-dh*x)
        vac=W*np.exp(-tot/T)
        fwd=vac*(1.-a)*(1.-b)
        inv=W*a*b
        pole=(abs(M[sl,None]-V)<5*width)
        for j,ev in enumerate((vac,fwd,inv,inv*pole)):
            accum[j]+=(ev.sum(),(ev*tot/2).sum())
    q0,_,_,N0,_=source.event_moments(T)
    reference=np.array([2*N0.sum(),np.sum(q0*N0)])
    assert np.max(abs(accum[0]/reference-1.))<2e-12
    return accum


def run(table,paths,crop_last):
    archive=[(Path(p).name,json.loads(Path(p).read_text())["snapshots"][-1]) for p in paths]
    if crop_last:
        name,s=archive[-1]
        archive.append((name+f": same-state crop q/T_R={crop_last:g}",cropped(s,crop_last)))
    results=[]
    for name,s in archive:
        T=s["T_visible_GeV"]
        H=np.sqrt(TabulatedRadiation(table,T).state(1.)[1]/3)/MPL
        n=s["represented_species_count_a3"][1]/s["scale_factor"]**3
        rr=[]
        pmax=s["lab_momentum_edges_GeV"][-1]
        lb=2*(np.hypot(pmax,m)-m)/T
        for order,polar,tail in (((8,20,16),6,max(55.,1.2*lb)),((12,35,22),9,max(100.,1.5*lb))):
            a=rates(T,s,order,polar,tail)
            vac,forward,rev,pole=a[:,0]
            assert 0.<=pole<=rev and 0.<=forward<=vac*(1.+2e-12)
            rr.append(dict(order=order,polar_order=polar,boost_excess_cut_over_T=tail,
                           direct_visible_vacuum_gross_number_rate_per_n_H=2*vac/(n*H),
                           direct_visible_forward_gross_number_rate_per_n_H=2*forward/(n*H),
                           direct_visible_reverse_gross_number_loss_per_n_H=2*rev/(n*H),
                           reverse_over_occupancy_corrected_forward=rev/forward,
                           fraction_of_reverse_events_within_five_mediator_widths=pole/rev if rev>0. else None))
        results.append(dict(snapshot=name,T_GeV=T,mediator_mass_GeV=V,
          maximum_represented_pair_invariant_mass_GeV=float(2*np.hypot(pmax,m)),rates=rr))
    return results


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("paths",nargs="+")
    p.add_argument("--crop-last",type=float)
    a=p.parse_args()
    print(json.dumps(dict(rows=run(a.table,a.paths,a.crop_last),
        scope="Gross direct stable pair processes into/from visible e, mu, tau; no heavy-parent channels, phase-space or cosmology evolution, or exterior particles; finite-width mediator and visible FD occupations. Rates are instantaneous; do not infer cumulative depletion."),
        indent=2,allow_nan=False))
