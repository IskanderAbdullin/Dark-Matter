"""Finite-T free-quark trace in a restricted hard timelike portal integrand."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import kve,roots_genlaguerre,roots_legendre

from conditional_relic_normalization_reference import conversion
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp
from timelike_portal_onsager_coverage_reference import factor,V,width

XEDGES=[4.,4.5,5.,6.,7.,8.,9.,10.,12.]
UEDGES=[.0125,.06,.15,.35,.65,1.,1.5,2.2,3.1,4.4,6.5,9.,13.,18.,24.,32.]


def gauss(edges,n):
    r,w=roots_legendre(n)
    e=np.asarray(edges)
    lo,hi=e[:-1,None],e[1:,None]
    return ((lo+hi+(hi-lo)*r)/2).ravel(),((hi-lo)*w/2).ravel()


def thermal_trace_over_vacuum(x,u):
    p=np.sqrt(u*(2*x+u))
    e=x+u
    # k_minus=(omega-k)/(2T). Rational evaluation avoids cancellation.
    minus=x*x/(2*(e+p))
    plus=(e+p)/2
    return 1-2*(np.logaddexp(0.,-minus)-np.logaddexp(0.,-plus))/p


def identities():
    x=np.array([4.,5.,7.,10.,20.])[:,None]
    u=np.array([.02,.4,1.,4.,10.,25.])[None,:]
    y,w=roots_legendre(36)
    q=np.sqrt(u*(2*x+u))
    e=x+u
    A=(e[...,None]+q[...,None]*y)/2
    B=(e[...,None]-q[...,None]*y)/2
    # f_F(A) f_F(B) divided by exp(-omega/T).
    fd=(1/(1+np.exp(-A))/(1+np.exp(-B)))@(w/2)
    kms=thermal_trace_over_vacuum(x,u)/(-np.expm1(-e))
    d1=float(np.max(np.abs(fd/kms-1.)))
    y,w=roots_genlaguerre(32,.5)
    J=np.sqrt(2*x+y)
    N=(J@w)/x[:,0]
    U=(((x+y)*J)@w)/x[:,0]
    exact=[kve(1,x[:,0]),2*kve(1,x[:,0])+x[:,0]*kve(0,x[:,0])]
    d2=float(np.max(np.abs(np.stack([N,U])/exact-1.)))
    assert max(d1,d2)<1e-11
    return dict(KMS_vs_angular_FD_relative_defect=d1,
                full_boost_MB_vs_closed_bessel_moments_relative_defect=d2)


def kernel(n,xedges=None):
    x,dx=gauss(XEDGES if xedges is None else xedges,n)
    u,du=gauss(UEDGES,n)
    x,dx=x[:,None],dx[:,None]
    q=np.sqrt(u*(2*x+u))
    e=x+u
    # Nc sum_{u,d,s} Q_f^2 = 2, with strictly massless quarks. In the
    # second convention, exp(-omega/T) replaces the full KMS*trace factor.
    R=2*np.stack([thermal_trace_over_vacuum(x,u),-np.expm1(-e)])
    # Dimensionless ds and boost measure: ds=T^2 2x dx;
    # integral (k/sqrt(s)) d omega = T integral (q/x) d u.
    W=2*x*dx*du*q/x/np.expm1(e)
    assert np.all(R>0.)
    return x[:,0],np.sum(R*W,axis=2),np.sum(R*(W*e),axis=2)


def local(T,kin):
    x,N,U=kin
    E=x*T
    s=E*E
    K=factor*E**3/((s-V*V)**2+(V*width)**2)
    ans=[]
    for mass in (m,M):
        b=K*(s+2*mass**2)*np.sqrt(np.maximum(1-4*mass**2/s,0.))*T**3
        ans.append(np.stack([2*N@b,T*U@b],axis=1))
    return np.stack(ans,axis=1)


def end(bg,normalization,kin,lo,hi,n):
    def densities(u):
        return (2*np.pi**2/45*np.exp(bg.entropy_log(u)),
                np.pi**2/30*np.exp(4*u+bg.rg(u)))
    sf,rf=densities(np.log(.001))
    edges=np.log([lo,hi])
    edges=np.unique(np.clip(np.r_[edges,bg.sg.x],*edges))
    us,ws=gauss(edges,n)
    z=0.
    for u,w in zip(us,ws):
        T=np.exp(u)
        entropy,rho=densities(u)
        dt=w*(3+bg.sg.derivative()(u))/(3*np.sqrt(rho/3)/MPL)
        z+=local(T,kin)*np.array([dt/entropy,dt*(sf/entropy)**(4/3)/rf])
    count=mp*normalization*(z[:,0,0]+2*z[:,1,0])
    energy=z[:,:,1].sum(axis=1)
    return dict(temperature_range_GeV=[lo,hi],
                count_omega_if_all_convert=count.tolist(),
                energy_floor_over_radiation_at_1MeV=energy.tolist())


def run(path):
    bg=TabulatedRadiation(path,.05)
    f=conversion(path)["Omega_h2_per_GeV_times_yield"]
    spans=[(.25,.34),(.25,.34862),(.25,.42),(.30,.34862),(.30,.42)]
    test=identities()
    rows=[]
    for n in (2,4,8):
        kin=kernel(n)
        rows.append(dict(nodes_per_panel=n,rows=[end(bg,f,kin,a,b,max(n,4)) for a,b in spans]))
    ds=[]
    for a,b in zip(rows[1]["rows"],rows[2]["rows"]):
        for k in ("count_omega_if_all_convert","energy_floor_over_radiation_at_1MeV"):
            ds.append(np.max(abs(np.array(a[k])/b[k]-1.)))
    assert max(ds)<1e-5
    return dict(case_order=["free three-flavour FD/KMS Born","matched free vacuum ratio with MB"],
        x_M_over_T_bounds=[XEDGES[0],XEDGES[-1]],
        u_omega_minus_M_over_T_bounds=[UEDGES[0],UEDGES[-1]],
        primary_trace_equation="https://arxiv.org/html/1910.09567#S2.E4",
        radiation_table_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        identities=test, max_last_refinement_defect=max(ds),controls=rows,
        scope="Restricted massless free-quark hard timelike Born test on the prescribed radiation history. Do not add to overlapping measured spectral bands. No interacting QCD, light-cone resummation, nonzero quark masses, positive inverse source or self-consistent dark kinetic history.")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("radiation_table")
    a=p.parse_args()
    print(json.dumps(run(a.radiation_table),indent=2,allow_nan=False))
