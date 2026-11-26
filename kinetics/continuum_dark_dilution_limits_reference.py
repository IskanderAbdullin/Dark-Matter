"""Conditional collisionless and reversible-adiabatic dilution limits."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, zeta

from continuum_dark_normal_capacity_reference import critical_T, fugacity, moments, rule, from_row
from sm_dark_entropy_exchange_reference import StandardModel
from three_species_decay_onsager_reference import m,mp


def initial_shells(snapshot,order=96):
    p=np.asarray(snapshot["lab_momentum_edges_GeV"])
    x,w=rule(order)
    pj=(p[:-1,None]*(1-x)+p[1:,None]*(1+x))/2
    out=[]
    for name,mass,g in (("scalar",mp,1),("fermion",m,4)):
        a,b=np.asarray(snapshot["entropy_variable_"+name])
        K=lambda px:px*px/(np.hypot(px,mass)+mass)
        z=a[:,None]+b[:,None]*(K(pj)-K(p[:-1,None]))/(K(p[1:,None])-K(p[:-1,None]))
        f=np.exp(-z)/(-np.expm1(-z)) if g==1 else expit(-z)
        u=g/(4*np.pi**2)*pj*pj*f*np.diff(p)[:,None]*w
        out.append((mass,pj,u))
    return out


def diagnostic(rows,snapshot,eos,extent):
    cap=from_row(rows,eos,20)
    shell=initial_shells(snapshot)
    ni=np.array([u.sum() for _,_,u in shell])
    assert np.max(abs(ni/np.array(cap['species_counts_GeV3'])-1.))<5e-12
    def kin(b):
        return sum(np.sum(u*(q/b)**2/(np.hypot(q/b,mass)+mass)) for mass,q,u in shell)/b**3
    Ki=kin(1.)
    assert abs(Ki/cap["retained_dark_kinetic_energy_GeV4"]-1.)<5e-12
    entropy0=cap["entropy_density_GeV3"]
    def at(b, emit=True):
        n=ni/b**3
        Tc=critical_T(n[0],20)
        ef=fugacity(Tc,m,n[1],False,20)
        Kc=moments(Tc,mp,0.,True)[1]+moments(Tc,m,ef,False,20)[1]
        Sc=moments(Tc,mp,0.,True,entropy=True)+moments(Tc,m,ef,False,entropy=True)
        S=entropy0/b**3
        if not emit: return S/Sc-1.
        def state(T, saturated=False):
            ep=0. if saturated else fugacity(T,mp,n[0],True,20)
            ef=fugacity(T,m,n[1],False,20)
            h=moments(T,mp,ep,True,entropy=True)+moments(T,m,ef,False,entropy=True)
            return h,float(np.exp(ep))
        if S>Sc*(1+1e-12):
            hi=2*Tc
            while state(hi)[0]<S: hi*=2
            Te=np.exp(brentq(lambda l:state(np.exp(l))[0]/S-1.,np.log(Tc*(1+1e-13)),np.log(hi)))
            hs,zp=state(Te)
            ground=0.
        elif S<Sc*(1-1e-12):
            lo=Tc*.5
            while state(lo,True)[0]>S:lo*=.5
            Te=np.exp(brentq(lambda l:state(np.exp(l),True)[0]/S-1.,np.log(lo),np.log(Tc)))
            hs,zp=state(Te,True)
            ground=max(0.,1.-moments(Te,mp,0.,True)[0]/n[0])
        else:
            Te,zp,hs,ground=Tc,1.,Sc,0.
        assert abs(hs/S-1)<1e-10
        r=kin(b)/Kc
        return dict(scale_factor_relative_to_snapshot=float(b),
                    collisionless_kinetic_energy_over_common_saturation=float(r),
                    reversible_adiabatic_entropy_over_common_saturation=float(S/Sc),
                    reversible_common_T_GeV=float(Te),
                    reversible_scalar_fugacity=zp,
                    reversible_equilibrium_condensate_fraction=float(ground),
                    saturation_T_GeV=float(Tc))

    out=[at(b) for b in extent]
    assert abs(out[0]["reversible_common_T_GeV"]/cap["common_T_GeV"]-1.)<2e-10
    flips=[(a,b) for a,b in zip(extent,extent[1:]) if at(a,False)*at(b,False)<0.]
    crossings=[at(np.exp(brentq(lambda l:at(np.exp(l),False),np.log(a),np.log(b))))
               for a,b in flips]
    # Independent very late all-nonrelativistic limit. The fermion gas is
    # Boltzmann at this sparse density/T ratio (test its fugacity).
    N=ni
    x=2.5*zeta(2.5)/zeta(1.5)
    eta=np.log(N[1]/(4*(m/mp)**1.5*N[0]/zeta(1.5)))
    assert np.exp(eta)<2e-6
    def late(z):
        return x*N[0]*z**1.5+N[1]*(2.5-eta+1.5*np.log(z))-entropy0
    if late(1.)>0:
        z=brentq(late,1e-5,1.)
        asymptote=1-z**1.5
        assert abs(out[-1]["reversible_equilibrium_condensate_fraction"]-asymptote)<2e-6
    else:
        asymptote=0.
    return dict(rows=out,crossings=crossings,
                nonrelativistic_asymptotic_condensate_fraction=float(asymptote),
                minimum_energy_ratio=min(x['collisionless_kinetic_energy_over_common_saturation'] for x in out),
                minimum_entropy_ratio=min(x['reversible_adiabatic_entropy_over_common_saturation'] for x in out))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("run")
    p.add_argument("snapshot")
    a=p.parse_args()
    row=json.loads(Path(a.run).read_text())["rows"][-1]
    snap=json.loads(Path(a.snapshot).read_text())["snapshots"][-1]
    assert row["nsteps"]==snap["steps"] and abs(row["final_visible_temperature_GeV"]-snap["T_visible_GeV"])<5e-12
    grid=np.r_[1.,np.geomspace(1.04,1e9,60)]
    print(json.dumps(diagnostic(row,snap,StandardModel(a.eos),grid),indent=2,allow_nan=False))
