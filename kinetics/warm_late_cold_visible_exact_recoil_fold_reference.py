"""Positive, endpoint-aligned FD/MB cold recoil fold for the direct vacuum channel."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from empty_lepton_pair_spectrum_reference import Source
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,Es,P,mp
from timelike_chemical_heat_eigenmodes_reference import leptons
from timelike_portal_onsager_coverage_reference import V,width,factor,gprod
from warm_late_cold_visible_offpole_dose_reference import log_tail,source
from warm_late_cold_visible_vacuum_dose_reference import run as base


DETAIL_MASS_EDGES=np.array([.2,.28,.4,.6,1.,np.inf])

def local(T, pc, n, detail=False):
    C,L,higher=limits(pc)
    y,w=roots_legendre(n)
    t=(y+1)/2
    tw=w/2
    # Resolve the narrow Boltzmann support above the exact total-energy gate.
    gate=L+M
    levels=np.array([0,.25,1,3,6,10,16,24,36,50,64,80])
    Q=(gate+T*(levels[:-1,None]+np.diff(levels)[:,None]*t)).ravel()
    dQ=(T*np.diff(levels)[:,None]*tw).ravel()
    r=np.zeros(2)
    by_flavor=np.zeros(len(leptons))
    by_mass=np.zeros(len(DETAIL_MASS_EDGES)-1)
    for q,dq in zip(Q,dQ,strict=True):
        if L>q/2:
            b=4*(L*(q-L)+M*M)
            disc=(b-4*M*q)*(b+4*M*q)
            assert disc>0
            vmax=(b+np.sqrt(disc))/2
            vmin=(2*M*q)**2/vmax
            low,high=np.sqrt([vmin,vmax])
        else:
            low,high=2*M,q
        assert low<high
        ends=np.unique(np.clip(np.r_[low,high,2*M,2*np.asarray(leptons),V,DETAIL_MASS_EDGES[np.isfinite(DETAIL_MASS_EDGES)] if detail else []],low,high))
        R=(ends[:-1,None]+np.diff(ends)[:,None]*t).ravel()
        dR=(np.diff(ends)[:,None]*tw).ravel()
        v=(R-2*M)*(R+2*M)/R**2
        pair=np.sqrt((q-R)*(q+R))
        half=pair*np.sqrt(v)/2
        Emin=q/2-half
        Emax=q/2+half
        a=np.maximum(L,Emin);b=np.minimum(L+higher,Emax)
        length=np.maximum(b-a,0.)
        Et=a[:,None]+length[:,None]*t
        u=(Et-q/2)/half[:,None]
        boost=np.sqrt((Et-M)*(Et+M))/M
        energy_low=Es/(Et/M+boost)+(mp**2/(Es+P))*boost
        decay=np.maximum(C-energy_low,0.)/(2*P*boost)
        assert np.all((decay>=0)&(decay<=1.))
        W=4*factor*(R**2+2*M*M)*R**3*np.sqrt(v)*pair*dR*dq/((R**2-V*V)**2+(V*width)**2)
        # Factor 2 for the two parent daughters is included in 4*factor.
        B0=np.zeros_like(R);B2=B0.copy();MB=B0.copy()
        W0=[]
        for ml in leptons:
            ok=R>2*ml
            fvel=(R[ok]-2*ml)*(R[ok]+2*ml)/R[ok]**2
            x=pair[ok,None]*np.sqrt(fvel)[:,None]*y/(2*T)
            e=q/(2*T)
            bl=1/(1+np.exp(-e+x)+np.exp(-e-x)+np.exp(-2*e))
            weight=(1+2*ml*ml/R[ok]**2)*np.sqrt(fvel)
            c2=(2/3)*v[ok]*fvel/(3-v[ok]-fvel+v[ok]*fvel/3)
            MB[ok]+=weight
            coeff0=weight*(bl@tw)
            coeff2=weight*c2*(bl@(tw*(3*y*y-1)/2))
            B0[ok]+=coeff0
            B2[ok]+=coeff2
            if detail:
                W0.append((ok,coeff0,coeff2))
        stat=B0[:,None]+B2[:,None]*(3*u*u-1)/2
        assert np.min(stat)>-1e-20 and np.max(stat-MB[:,None])<2e-14
        measure=W*np.exp(-q/T)*length/(2*half)
        d0=decay@tw
        d2=(decay*(3*u*u-1)/2)@tw
        fd=measure*((stat*decay)@tw)
        r+=np.array([fd.sum(),measure@(MB*d0)])
        if detail:
            for i,(ok,c0,c2) in enumerate(W0):
                by_flavor[i]+=measure[ok]@(c0*d0[ok]+c2*d2[ok])
            k=np.searchsorted(DETAIL_MASS_EDGES,R,side="right")-1
            by_mass+=np.bincount(k,weights=fd,minlength=len(by_mass))
    x=(gate+80*T)/T
    bw=(np.hypot(V*V,V*width)+V*V)/(2*(V*width)**2)
    c=2*12*3*gprod*gprod*bw/(8*np.pi)/(120*np.pi**4)
    tail=float(np.exp(np.log(c)+6*np.log(T)+log_tail(x)))
    assert 0<r[0]<=r[1] and tail<1e-9*r[1], (T, pc, r, tail)
    if detail:
        assert abs(by_flavor.sum()/r[0]-1)<3e-14
        assert abs(by_mass.sum()/r[0]-1)<3e-14
        return np.r_[r,tail,by_flavor,by_mass]
    return np.r_[r,tail]


def run(eos,archive):
    bg=Background(eos)
    y=base(eos,archive)["comparison_archive_held_equality_cold_yield"]
    rows=[]
    # Below 12 MeV, the published whole-bin upper is already negligible;
    # here give a positive numerical estimate of only the 12–TR band.
    for start in (.025,.035,.05):
        rad=TabulatedRadiation(eos,start)
        results=[]
        for n,nt in ((8,12),(12,20)):
            x,w=roots_legendre(nt)
            a,b=np.log([.012,start])
            L=(b-a)/2
            T=np.exp((a+b)/2+L*x)
            ys=[]
            for t in T:
                A=bg.scale(t)
                tj,rho=rad.state(A/bg.scale(start))
                assert abs(tj/t-1)<3e-12
                entropy=2*np.pi**2/45*np.exp(bg.gs(np.log(t)))*t**3
                dt_dlogT=(3+bg.gs.derivative()(np.log(t)))/3
                pc=.1*bg.T0/A
                rr=local(t,pc,n)
                # Analytic bin envelope is an independent sanity check.
                assert rr[1]<2*np.exp(source(t,pc,[0,1,2,4,7,10,np.inf]))
                ys.append(rr*dt_dlogT/(entropy*np.sqrt(rho/3)/MPL))
            Y=L*(w@np.asarray(ys))
            results.append(dict(kinetic_order=n,temperature_order=nt,FD_yield=float(Y[0]),MB_yield=float(Y[1]),
                                integrated_sample_of_local_tail_upper=float(Y[2]),FD_to_held=float(Y[0]/y),MB_to_held=float(Y[1]/y)))
        assert abs(results[0]["FD_yield"]/results[1]["FD_yield"]-1)<.03
        rows.append(dict(start_temperature_GeV=start,controls=results))
    # Recompute a distinct 50 MeV fold with explicit invariant bins.
    start=.05
    rad=TabulatedRadiation(eos,start)
    x,w=roots_legendre(20)
    a,b=np.log([.012,start]);h=(b-a)/2
    detail_y=[]
    for t in np.exp((a+b)/2+h*x):
        A=bg.scale(t)
        _,rho=rad.state(A/bg.scale(start))
        entropy=(2*np.pi**2/45)*np.exp(bg.gs(np.log(t)))*t**3
        dtime=(3+bg.gs.derivative()(np.log(t)))/3
        detail_y.append(local(t,.1*bg.T0/A,12,detail=True)*dtime/(entropy*np.sqrt(rho/3)/MPL))
    detail_sum=h*(w@np.array(detail_y))
    assert abs(detail_sum[0]/rows[-1]["controls"][-1]["FD_yield"]-1)<5e-4
    detailed=dict(flavors=["electron","muon","tau"], FD_fraction_by_incoming_flavor=(detail_sum[3:6]/detail_sum[0]).tolist(),
                  invariant_mass_edges_GeV=[float(x) if np.isfinite(x) else None for x in DETAIL_MASS_EDGES],
                  FD_fraction_by_pair_invariant_mass=(detail_sum[6:]/detail_sum[0]).tolist(),FD_total_yield=float(detail_sum[0]))
    checks=[]
    for T in (.035,.05):
        pc=.1*bg.T0/bg.scale(T)
        e=np.hypot(mp,pc)
        v=float(local(T,pc,12)[0])
        E,W=Source((24,64,32),tail=80.,dark_mass=M).particles(T,72)
        beta=np.sqrt((E-M)*(E+M))/M
        lo=Es/(E/M+beta)+(mp**2/(Es+P))*beta
        alt=float(W@(np.maximum(e-lo,0.)/(2*P*beta)))
        assert abs(alt/v-1)<4e-4
        checks.append(dict(T_GeV=T,total_energy_first=v,alternate_invariant_mass_first=alt))
    return dict(detailed_50_MeV=detailed,independent_local_chart_checks=checks,scope="Positive numerical FD and MB estimates, not certified full upper. Open unpolarized vacuum decay; "
                      "direct thermal charged leptons into heavy pair only; 12–TR MeV, adiabatic and no later scattering.",
                held_equality_cold_yield=y,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("archive",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
