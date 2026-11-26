"""Conditional common-temperature continuum FD/BE capacity test.

Conserved inputs are the separate represented scalar and stable-fermion
counts and their combined kinetic energy. A continuum common equilibrium is
a diagnostic, not a claim that the dynamical spectra reach it.
"""
import argparse
import json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre, zeta

from sm_dark_entropy_exchange_reference import StandardModel
from three_species_decay_onsager_reference import m, mp


@lru_cache(maxsize=8)
def rule(order):
    return roots_legendre(order)


def moments(T, mass, eta, bose, order=20, entropy=False):
    """Number and kinetic energy, eta=(mu-m)/T, integrated to negligible tails."""
    assert T>0 and mass>=0 and (not bose or eta<=0.)
    end=max(0.,eta)+72.
    edges=np.unique(np.clip(np.r_[0., .03, .12, 1., 4., 12., 24., 60.,
                              eta+np.array([-18.,-7.,0.,7.,18.,36.]), end],0.,end))
    lo,hi=np.sqrt(edges[:-1]),np.sqrt(edges[1:])
    x,w=rule(order)
    u=(lo[:,None]*(1-x)+hi[:,None]*(1+x))/2
    K=T*u*u
    E=mass+K
    p=np.sqrt(K*(K+2*mass))
    f=np.exp(-u*u+eta)/(-np.expm1(-u*u+eta)) if bose else expit(eta-u*u)
    g=1 if bose else 4
    measure=g/(2*np.pi**2)*p*E*(T*u*(hi-lo)[:,None]*w)
    if entropy:
        lam=u*u-eta
        hb=lam*f-np.log(-np.expm1(-lam)) if bose else lam*f+np.logaddexp(0.,-lam)
        return float(np.sum(measure*hb))
    return np.array([np.sum(measure*f), np.sum(measure*f*K)])


def fugacity(T, mass, density, bose, order):
    hi=0. if bose else 1.
    def fun(q):
        return moments(T,mass,q,bose,order)[0]/density-1.
    if bose and fun(hi)<0.:
        if fun(hi)>-1e-10: return hi
        raise ValueError("normal Bose state is above the saturation count")
    if not bose:
        while fun(hi)<0.:
            hi=2*hi+2
    lo=-32.
    while fun(lo)>0.: lo*=2
    return brentq(fun,lo,hi,xtol=1e-12)


def critical_T(n,order):
    typ=np.cbrt(n*np.pi**2/zeta(3.))
    return float(np.exp(brentq(lambda z:moments(np.exp(z),mp,0.,True,order)[0]/n-1.,
                              np.log(typ)-24,np.log(typ)+1e-10)))


def solve(n, kinetic, order=20):
    n=np.asarray(n,dtype=float)
    assert np.all(n>0.) and kinetic>0.
    def thermal(T,saturated):
        phi=0. if saturated else fugacity(T,mp,n[0],True,order)
        chi=fugacity(T,m,n[1],False,order)
        ns,ks=moments(T,mp,phi,True,order)
        nf,kf=moments(T,m,chi,False,order)
        return np.array([ns,nf]),ks+kf,np.array([phi,chi])
    Tc=critical_T(n[0],order)
    _,kc,eta=thermal(Tc,True)
    excess=(kinetic>=kc)
    if excess and kinetic/kc-1.>1e-12:
        tlo=Tc*(1+1e-13); thi=Tc*2.
        while thermal(thi,False)[1]<kinetic: thi*=2.
        T=np.exp(brentq(lambda l:thermal(np.exp(l),False)[1]/kinetic-1.,np.log(tlo),np.log(thi)))
    elif not excess and 1.-kinetic/kc>1e-12:
        tlo=Tc*.25
        # For this independently tested branch choose an input above the
        # fixed-number zero-temperature fermion kinetic energy.
        while thermal(tlo,True)[1]>kinetic:
            tlo*=.5
            if tlo<Tc*1e-5: raise ValueError("kinetic input below tested thermal bracket")
        T=np.exp(brentq(lambda l:thermal(np.exp(l),True)[1]/kinetic-1.,np.log(tlo),np.log(Tc)))
    else: T=Tc
    ns,k,eta=thermal(T,not excess or T==Tc)
    n0=max(0.,n[0]-ns[0]) if not excess else 0.
    assert max(abs(k/kinetic-1.),abs(ns[1]/n[1]-1.),abs((ns[0]+n0)/n[0]-1.))<2e-8
    rest=mp*n[0]+m*n[1]
    return dict(saturation_T_GeV=float(Tc),common_T_GeV=float(T),normal_common_phase=bool(excess),
                kinetic_energy_over_saturation=float(kinetic/kc),
                total_energy_over_saturation=float((kinetic+rest)/(kc+rest)),
                scalar_kinetic_fugacity=float(np.exp(eta[0])),
                fermion_kinetic_fugacity=float(np.exp(eta[1])),
                entropy_density_GeV3=float(moments(T,mp,eta[0],True,order,entropy=True)
                                          +moments(T,m,eta[1],False,order,entropy=True)),
                continuum_common_state_condensate_fraction=float(n0/n[0]),
                closure_kinetic_relative_defect=float(abs(k/kinetic-1.)))


def from_row(row,eos,order):
    a=row['final_a']
    s=eos.state(row['initial_visible_T_GeV'])[2]
    count=np.array(row['final_counts_div_initial_visible_entropy'])*s/a**3
    # Reconstruct the dimensional retained dark energy from two independently
    # archived ledgers: total source less both tails, and final E/input.
    full=np.array(row['comoving_source_channel_count_and_total_energy'])
    tails=np.array(row['comoving_source_channel_count_and_total_energy_tails'])
    keep=(full-tails.sum(axis=2))[:,1].sum()+sum(row['first_pulse_energy_a4'])
    E=keep*sum(row['final_energy_over_scale'])/a**4
    kinetic=E-(np.array([mp,m])@count)
    return dict(nsteps=row["nsteps"],actual_visible_T_GeV=row["final_visible_temperature_GeV"],
                species_counts_GeV3=count.tolist(),retained_dark_kinetic_energy_GeV4=float(kinetic),
                **solve(count,kinetic,order))


def shape_metrics(cap, snapshot, order=96):
    """Compare archived kinetic shells with the continuum common normal state."""
    assert cap["normal_common_phase"]
    p=np.asarray(snapshot["lab_momentum_edges_GeV"])
    T=cap["common_T_GeV"]
    x,w=rule(order)
    q=(p[:-1,None]*(1-x)+p[1:,None]*(1+x))/2
    width=np.diff(p)[:,None]/2
    hsum=0.
    output=[]
    for j,(mass,bose,fug,name) in enumerate(((mp,True,cap["scalar_kinetic_fugacity"],"scalar"),
                                             (m,False,cap["fermion_kinetic_fugacity"],"fermion"))):
        g=1 if bose else 4
        a,b=np.array(snapshot["entropy_variable_"+name])
        kinetic=lambda s:s*s/(np.hypot(s,mass)+mass)
        low,high=kinetic(p[:-1]),kinetic(p[1:])
        k=kinetic(q)
        z=a[:,None]+b[:,None]*(k-low[:,None])/(high-low)[:,None]
        r=k/T-np.log(fug)
        f=np.exp(-z)/(-np.expm1(-z)) if bose else expit(-z)
        ref=np.exp(-r)/(-np.expm1(-r)) if bose else expit(-r)
        ss=z*f-np.log(-np.expm1(-z)) if bose else z*f+np.logaddexp(0.,-z)
        u=g/(2*np.pi**2)*q*q*width*w
        n=cap["species_counts_GeV3"][j]
        N=float(np.sum(u*f))
        Nref=float(np.sum(u*ref))
        assert abs(N/n-1.)<5e-12 and Nref<=n*(1+3e-12)
        outside=max(0.,1-Nref/n)
        metric=.5*(np.sum(u*abs(f-ref))/n+outside)
        hsum+=np.sum(u*ss)
        output.append(dict(species=name,closed_count_defect=abs(N/n-1.),
                           continuum_common_count_outside_original_mesh_fraction=outside,
                           number_normalized_total_variation=float(metric)))
    assert hsum<=cap["entropy_density_GeV3"]*(1+2e-12)
    return dict(actual_entropy_over_continuum_maximum=float(hsum/cap["entropy_density_GeV3"]), species=output)

def selftest():
    T=.07
    light=moments(T,0.,0.,True)
    n_theory=T**3*zeta(3)/np.pi**2
    k_theory=T**4*np.pi**2/30
    massless_err=max(abs(light/np.array([n_theory,k_theory])-1.))
    T=1e-6; mass=.2
    heavy=moments(T,mass,0.,True)
    B=(mass*T/(2*np.pi))**1.5
    nonrel_err=max(abs(heavy/(B*np.array([zeta(1.5),1.5*T*zeta(2.5)]))-1.))
    assert massless_err<1e-12 and nonrel_err<2e-5
    T=.01; Np,Kp=moments(T,mp,-.45,True);Nf,Kf=moments(T,m,-1.5,False)
    normal=solve([Np,Nf],Kp+Kf)
    assert normal['normal_common_phase'] and abs(normal['common_T_GeV']/T-1.)<1e-8
    Np,Kp=moments(T,mp,0.,True)
    test=solve([1.25*Np,Nf],Kp+Kf)
    assert not test['normal_common_phase'] and abs(test['continuum_common_state_condensate_fraction']-.2)<1e-8
    return dict(massless_exact_relative_defect=float(massless_err),
                nonrelativistic_asymptotic_relative_defect=float(nonrel_err),
                manufactured_normal=normal,manufactured_condensed=test)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("runs",nargs="*")
    p.add_argument("--snapshots",type=str,nargs="+")
    a=p.parse_args()
    eos=StandardModel(a.eos)
    rows=[]
    if a.snapshots is not None: assert len(a.snapshots)==len(a.runs)
    for idx,f in enumerate(a.runs):
        row=json.loads(Path(f).read_text())["rows"][-1]
        measured=from_row(row,eos,20)
        fine=from_row(row,eos,32)
        defect=max(abs(measured[k]/fine[k]-1.) for k in ("common_T_GeV","kinetic_energy_over_saturation","scalar_kinetic_fugacity"))
        assert defect<1e-8
        if a.snapshots is not None:
            sn=json.loads(Path(a.snapshots[idx]).read_text())["snapshots"][-1]
            assert sn["steps"]==row["nsteps"]
            assert abs(sn["T_visible_GeV"]-row["final_visible_temperature_GeV"])<5e-12
            measured["kinetic_shape_comparison"]=shape_metrics(measured,sn)
        rows.append(dict(run=Path(f).name,quadrature_refinement_relative_defect=defect,**measured))
    print(json.dumps(dict(controls=selftest(),rows=rows),indent=2,allow_nan=False))
