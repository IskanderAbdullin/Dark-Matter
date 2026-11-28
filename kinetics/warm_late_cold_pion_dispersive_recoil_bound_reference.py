"""Transfer positive vacuum dipion moments into the cold recoil yield."""
import argparse,json
from pathlib import Path

import numpy as np
from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import Es,M,P,mp
from timelike_portal_onsager_coverage_reference import V,width,mpi,factor
from warm_babar_dispersive_gap_reference import POINT,mu_kernel_lower
from warm_late_cold_visible_recoil_box_upper_reference import log_slope_minimum
from warm_light_occupation_stage_audit_reference import read


def per_sigma(T,pc,l,h,nq,nx):
    """Upper bound on d(weighted rate)/ds per unit sigma in a whole s bin."""
    C,L,_=limits(pc)
    S0=L+M
    y=np.r_[np.linspace(0,12,nq+1),np.linspace(12,36,nq//2+1)[1:],np.linspace(36,80,nq//2+1)[1:]]
    Q=S0+T*y
    Ql,Qh=Q[:-1,None],Q[1:,None]
    lo,hi=np.sqrt(l)[None,:],np.sqrt(h)[None,:]
    up=np.sqrt(np.maximum(Qh-lo,0)*(Qh+lo))
    dn=np.sqrt(np.maximum(Ql-hi,0)*(Ql+hi))
    v0=(l-4*M*M)/l;v1=(h-4*M*M)/h
    D0=dn*np.sqrt(v0)/2
    D1=up*np.sqrt(v1)/2
    x=np.linspace(-1,1,nx+1)[1:]
    Eh=np.minimum(Qh[:,:,None]/2+np.where(x>=0,D1[:,:,None],D0[:,:,None])*x,(Qh-M)[:,:,None])
    gam=np.maximum(L,np.minimum(Eh,M*Es/C))/M
    z=np.sqrt((gam-1)*(gam+1))
    el=Es/(gam+z)+(mp**2/(Es+P))*z
    A=np.where(Eh>L,np.maximum(C-el,0.)/(2*P*z),0.).mean(axis=2)
    # Vacuum unpolarized MB distribution is isotropic in incident CoM.
    # Bound the B-E enhancement before any nonnegative angular weights.
    E1=np.maximum(mpi,Ql/2-up*np.sqrt(np.maximum(0,1-4*mpi*mpi/h))/2)
    E2=np.maximum(E1,Ql-E1)
    BE=1./((-np.expm1(-E1/T))*(-np.expm1(-E2/T)))
    assert np.max(BE)<=(-1/np.expm1(-mpi/T))**2*(1+1e-14)
    b=(np.clip(V*V,l,h)-V*V)**2+(V*width)**2
    X=2*factor*(h+2*M*M)*h*h*np.sqrt(v1)/(POINT*b)
    qi=T*np.exp(-Ql/T)*(-np.expm1(-(Qh-Ql)/T))
    common=(A*up*qi*BE).sum(axis=0)
    tail=T*(Q[-1]+T)*np.exp(-Q[-1]/T)*(-1/np.expm1(-mpi/T))**2
    return X*(common+tail)


def phase_bound(bg,lo,hi,nt,nq,nx,*,qcut_today=None):
    # Whole time interval upper, split 12-50 more tightly. The final
    # coarse endpoint includes all the range to 0.1 MeV.
    T=np.r_[np.geomspace(.05,.012,nt+1),.0001]
    A=np.array([bg.scale(t) for t in T])
    z=np.log(T[1:])
    entropy=(2*np.pi**2/45)*np.exp(bg.gs(z))*T[1:]**3
    H=np.sqrt(np.pi**2/90*np.exp(bg.gr(z)))*T[1:]**2/MPL
    w=np.diff(np.log(A))/(entropy*H)
    if qcut_today is None:
        qcut_today=.1*bg.T0
    out=np.asarray([per_sigma(t,qcut_today/a,lo,hi,nq,nx) for t,a in zip(T[:-1],A[:-1],strict=True)])
    return w@out


def run(eos,babar,theory,held):
    bg=Background(eos)
    certify=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    K=read(held)["held_equality_cold_yield"]
    b=read(babar)
    t=read(theory)
    B=next(r["input_whole_dispersion_moment"] for r in b["moment_variants"]
           if r["data"]=="BABAR 2026" and r["linear_error_multiple"]==2)
    C=t["variants"][-1]["theory_full_positive_moment"]
    out=[]
    for name,l,h,cap in (
      ("lower threshold theory",float(t["threshold_s_GeV2"][0]),.09,C),
      ("BABAR low gap",.09,.10,B)):
        rows=[]
        for ns,nt,nq,nx in ((12,24,20,18),(24,48,36,30),(32,96,48,42)):
            edge=np.linspace(l,h,ns+1)
            weights=phase_bound(bg,edge[:-1],edge[1:],nt,nq,nx)
            mu=mu_kernel_lower(edge[1:])/(4*np.pi**3)
            sup=weights/mu
            imax=int(np.argmax(sup))
            obs=float(sup[imax]*cap)
            rows.append(dict(s_bins=ns,T_bins=nt,Q_base_bins=nq,parent_angle_bins=nx,
               max_yield_per_dispersion_unit=float(sup[imax]),cold_yield_upper=obs,held_fraction_upper=obs/K,
               controlling_s_GeV2=edge[imax:imax+2].tolist()))
        assert rows[-1]["cold_yield_upper"]<rows[0]["cold_yield_upper"]
        out.append(dict(name=name,s_GeV2=[l,h],assumed_positive_dispersion_cap=cap,controls=rows))
    return dict(scope="Conditional positive BABAR/theory caps on separate unmatched vacuum charged-pion intervals. "
             "Adiabatic 50 to .1 MeV, exact recoil inequality, bounded physical massive pion Bose occupation; "
             "not a statistical combination, unmeasured high bands or full thermal hadronic source.",
             entropy_radiation_log_slope_lower=certify,held_equality_cold_yield=K,bands=out,
             sum_of_conditional_outside_KLOE_cold_yield_upper=sum(b["controls"][-1]["cold_yield_upper"] for b in out))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","babar","theory","held"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.babar,a.theory,a.held),indent=2,allow_nan=False))
