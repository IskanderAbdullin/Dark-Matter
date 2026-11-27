"""Enclose the direct FD source using total-energy, invariant-energy and lab-angle boxes."""

import argparse
import json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import Es,M,P,mp
from timelike_chemical_heat_eigenmodes_reference import leptons
from timelike_portal_onsager_coverage_reference import V,width,factor,gprod
from warm_late_cold_visible_offpole_dose_reference import log_tail
from warm_late_cold_visible_exact_recoil_fold_reference import local
from warm_late_cold_visible_vacuum_dose_reference import run as base


def log_slope_minimum(interp,offset):
    a,b,c,_=interp.c
    h=np.diff(interp.x)
    x=np.clip(np.divide(-b,3*a,where=a!=0,out=np.zeros_like(a)),0,h)
    value=lambda d: offset+3*a*d*d+2*b*d+c
    out=float(min(value(0).min(),value(h).min(),value(x).min()))
    assert out>0
    return out


def rate(T,cut,nq,nr,nx):
    C,L,span=limits(cut)
    Q0=L+M
    y=np.r_[np.linspace(0,12,nq+1),np.linspace(12,36,nq//2+1)[1:],np.linspace(36,80,nq//2+1)[1:]]
    Q=Q0+T*y
    # Choose fixed boxes in each total-energy bin; the upper R wall
    # can conservatively include a triangle outside the actual domain.
    j=np.linspace(0,1,nr+1)
    R=2*M+(Q[1:,None]-2*M)*j
    rl,rh=R[:,:-1],R[:,1:]
    Ql,Qh=Q[:-1,None],Q[1:,None]
    pair_hi=np.sqrt((Qh-rl)*(Qh+rl))
    pair_lo=np.sqrt(np.maximum(Ql-rh,0)*(Ql+rh))
    vlo=(rl-2*M)*(rl+2*M)/rl**2
    vhi=(rh-2*M)*(rh+2*M)/rh**2
    delta_hi=pair_hi*np.sqrt(vhi)/2
    delta_lo=pair_lo*np.sqrt(vlo)/2
    x=np.linspace(-1,1,nx+1)[1:]
    emax=np.minimum(Qh[:,:,None]/2+np.where(x>=0,delta_hi[:,:,None],delta_lo[:,:,None])*x,
                    (Qh-M)[:,:,None])
    s=np.maximum(L,np.minimum(emax,Es*M/C))/M
    b=np.sqrt((s-1)*(s+1))
    elow=Es/(s+b)+(mp*mp/(Es+P))*b
    angular=np.where(emax>L,np.maximum(C-elow,0.)/(2*P*b),0.).mean(axis=2)
    xm=np.linspace(-1,1,nx+1)[:-1]
    emin=Ql[:,:,None]/2+np.where(xm>=0,delta_lo[:,:,None],delta_hi[:,:,None])*xm
    def exact_decay_at(ee):
        g=np.maximum(L,ee)/M
        s=np.sqrt((g-1)*(g+1))
        e=Es/(g+s)+(mp*mp/(Es+P))*s
        return np.where((ee>L)&(ee<L+span),np.maximum(C-e,0.)/(2*P*s),0.)
    angle_lower=np.minimum(exact_decay_at(emin),exact_decay_at(emax)).mean(axis=2)
    assert np.all(angle_lower<=angular+2e-14)
    ml=np.asarray(leptons)[:,None,None]
    lep=((1+2*(ml/rh)**2)*np.sqrt(np.maximum(1-4*(ml/rh)**2,0.))).sum(axis=0)
    # In the incident lepton pair the sum is >= Ql and each energy
    # is >= the lab lower endpoint below. The convex exponential
    # makes the FD denominator largest at the most unequal permissible
    # split; use that inequality *before* angular averaging.
    ell=np.maximum(ml,Ql/2-pair_hi*np.sqrt(np.maximum(1-4*(ml/rh)**2,0.))/2)
    other=np.maximum(ell,Ql-ell)
    factor_FD=(1+np.exp(-ell/T)+np.exp(-other/T)+np.exp(-np.maximum(Ql,2*ell)/T))
    assert np.all(factor_FD<=(1+np.exp(-ml/T))**2+2e-14)
    lep_low=(((1+2*(ml/rl)**2)*np.sqrt(np.maximum(1-4*(ml/rl)**2,0.)))/factor_FD).sum(axis=0)
    bw_far=np.maximum((rl**2-V*V)**2,(rh**2-V*V)**2)+(V*width)**2
    bw_den=(np.clip(V*V,rl*rl,rh*rh)-V*V)**2+(V*width)**2
    qint=T*np.exp(-Ql/T)*(-np.expm1(-(Qh-Ql)/T))
    upper=4*factor*(rh**2+2*M*M)*rh**3*np.sqrt(vhi)*pair_hi/bw_den
    wrapped=float(np.sum(qint*(rh-rl)*upper*lep*angular))
    lower=4*factor*(rl**2+2*M*M)*rl**3*np.sqrt(vlo)*pair_lo/bw_far
    wrapped_low=float(np.sum(qint*(rh-rl)*lower*lep_low*angle_lower*(rh<=Ql)))
    assert 0<=wrapped_low<wrapped
    sb=(np.hypot(V*V,V*width)+V*V)/(2*(V*width)**2)
    const=2*12*3*gprod*gprod*sb/(8*np.pi)/(120*np.pi**4)
    remainder=float(np.exp(np.log(const)+6*np.log(T)+log_tail(Q[-1]/T)))
    assert wrapped>0 and remainder<1e-9*wrapped
    return np.array([wrapped,remainder,wrapped_low])


def run(eos,archive):
    bg=Background(eos)
    scales=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    before=base(eos,archive)
    Y=before["comparison_archive_held_equality_cold_yield"]
    late=10**before["rows"][2]["refinements"][-1]["log10_vacuum_soft_yield_upper"]
    local_controls=[]
    for t in (.025,.035,.05):
        pc=.1*bg.T0/bg.scale(t)
        fold=local(t,pc,12)
        u=rate(t,pc,48,64,36)
        assert u[0]+u[1]>=fold[1]>=fold[0]>=u[2]
        local_controls.append(dict(temperature_GeV=t,exact_sampled_local_MB=float(fold[1]),
                                  box_local_upper=float(u[0]+u[1]), box_local_FD_lower=float(u[2])))
    results=[]
    for top in (.025,.035,.05):
        rows=[]
        for nt,nq,nr,nx in ((24,16,20,12),(48,32,40,24),(96,48,64,36),(192,72,96,54)):
            T=np.geomspace(top,.012,nt+1)
            A=np.array([bg.scale(u) for u in T])
            s=(2*np.pi**2/45)*np.exp(bg.gs(np.log(T[1:])))*T[1:]**3
            H=np.sqrt((np.pi**2/90)*np.exp(bg.gr(np.log(T[1:]))))*T[1:]**2/MPL
            # Occupation MB bound grows with earlier T at fixed phases;
            # larger injection cut only adds positive decays.
            q=np.array([rate(t,.1*bg.T0/a,nq,nr,nx) for t,a in zip(T[:-1],A[:-1],strict=True)])
            u=np.diff(np.log(A))@(q[:,:2]/(s*H)[:,None])
            v=np.array([rate(t,.1*bg.T0/a,nq,nr,nx)[2] for t,a in zip(T[1:],A[1:],strict=True)])
            sht=(2*np.pi**2/45)*np.exp(bg.gs(np.log(T[:-1])))*T[:-1]**3
            Hht=np.sqrt(np.pi**2/90*np.exp(bg.gr(np.log(T[:-1]))))*T[:-1]**2/MPL
            lower=float(np.diff(np.log(A))@(v/(sht*Hht)))
            full=float(u.sum()+late)
            assert 0<=lower<full
            rows.append(dict(temperature_bins=nt,total_energy_base_bins=nq,invariant_bins=nr,parent_angle_bins=nx,
                upper_12_to_start=float(u[0]),truncated_tail_upper=float(u[1]),early_plus_low_T_check_upper=full,
                upper_fraction_of_held_yield=full/Y, lower_FD_yield=lower, lower_fraction_of_held_yield=lower/Y))
        assert rows[-1]["upper_fraction_of_held_yield"]<rows[0]["upper_fraction_of_held_yield"]
        results.append(dict(start_T_GeV=top,controls=rows))
    return dict(local_controls=local_controls,scope="Analytical interval ceilings and floors for thermal charged-lepton heavy-pair production and open unpolarized "
                      "vacuum decay, monotone EOS photon-time background and no scatter-in. Ending at .1 MeV.",
                log_entropy_and_radiation_slope_lower=scales,held_yield=Y,reused_sub12_upper=late,rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("archive",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
