"""Two distinct Dirac chi masses scatter on ONE evolving finite lepton bath.

Shared-grid on-shell finite-vector 2<->2 Pauli events. The identical
Fermi-logit interpolation/deposition proves Sdot>=0 and exact
Nchi1/Nchi2/Nlepton and total energy. Only synthetic small-grid tests;
NO real phi Bose, finite SM weak plasma, physical gamma_el, FRW, or relic.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit,logit,xlogy
from numpy.polynomial.legendre import leggauss
from kinetics.independent_chi_lepton_5d_quadrature import radial_nodes


def cell_weights(grid,mass,g=4):
    d=grid[1]-grid[0]
    left=np.maximum(mass,grid-d/2)
    right=np.minimum(grid[-1],grid+d/2)
    return g/(6*math.pi**2)*(np.maximum(right*right-mass*mass,0)**1.5-
                                np.maximum(left*left-mass*mass,0)**1.5)


def interp_fermi(grid,f):
    odds=logit(f)
    def func(E):
        z=np.asarray(E)
        if np.any(z<grid[0]-1e-10) or np.any(z>grid[-1]+1e-10):
            raise ValueError("outside physical shared energy mesh")
        return expit(np.interp(z,grid,odds))
    return func


def deposit(arr,grid,energy,amount):
    vals=np.atleast_1d(np.asarray(energy,float)).ravel()
    dq=np.atleast_1d(np.asarray(amount,float)).ravel()
    if vals.shape!=dq.shape:raise ValueError("source dimension mismatch")
    step=grid[1]-grid[0]
    if np.min(vals)<grid[0]-1e-10 or np.max(vals)>grid[-1]+1e-10:
        raise ValueError("outgoing energy above shared cutoff")
    z=np.clip((vals-grid[0])/step,0,len(grid)-1-1e-13)
    i=np.floor(z).astype(int)
    frac=z-i
    np.add.at(arr,i,(1-frac)*dq)
    np.add.at(arr,i+1,frac*dq)


def system(masses=(.1,.02),mlepton=0.,bins=12,cut=12.):
    if not(len(masses)==2 and 0<min(masses) and max(masses)>min(masses)
           and mlepton>=0 and bins>=4 and cut>0):raise ValueError
    upper=max(masses)+mlepton+2*cut
    chi_grid=[np.linspace(M,upper,bins+1) for M in masses]
    lg=np.linspace(mlepton,upper,bins+1)
    # Distinct synthetic initial spectra; NOT χ2 produced by a real phi.
    initial_chi=[expit(-chi_grid[0]/.9+.15*np.exp(-chi_grid[0]/.9)),
                 expit(-chi_grid[1]/.8-.10*np.exp(-chi_grid[1]/.5))]
    initial_lepton=expit(-lg+.05*np.exp(-lg))
    return {"masses":tuple(masses),"mlepton":mlepton,"bins":bins,"cut":cut,
            "chi_grid":chi_grid,"lg":lg,
            "chi_w":[cell_weights(gr,M) for gr,M in zip(chi_grid,masses)],
            "lw":cell_weights(lg,mlepton),
            "initial_chi":initial_chi,"initial_lepton":initial_lepton}


def ledger(cfg,chi,lep):
    Ns=[float(np.dot(w,f)) for w,f in zip(cfg["chi_w"],chi)]
    Nl=float(np.dot(cfg["lw"],lep))
    H=sum(float(-np.dot(w,xlogy(f,f)+xlogy(1-f,1-f)))
          for w,f in zip(cfg["chi_w"],chi))
    H+=float(-np.dot(cfg["lw"],xlogy(lep,lep)+xlogy(1-lep,1-lep)))
    E=sum(float(np.dot(w,g*f)) for w,g,f in zip(cfg["chi_w"],cfg["chi_grid"],chi))
    E+=float(np.dot(cfg["lw"],cfg["lg"]*lep))
    return {"Nchi1":Ns[0],"Nchi2":Ns[1],"Nlepton":Nl,
            "Etotal":E,"entropy":H}


def rhs(cfg,chi,lep,*,nr=4,nmu=4,nc=8,nphi=8,G=.001,mv=.1,
        radial_scheme='uniform_p',segments=1,incoming_cut=None):
    if len(chi)!=2 or any(np.any((f<=0)|(f>=1)) for f in [*chi,lep]):
        raise ValueError("strict Fermi interior required")
    quadrature_cut=cfg['cut'] if incoming_cut is None else incoming_cut
    if not (math.isfinite(quadrature_cut) and 0<quadrature_cut<=cfg['cut']):
        raise ValueError('incoming physical cutoff must not exceed energy mesh coverage')
    lepf=interp_fermi(cfg["lg"],lep)
    dchi=[np.zeros_like(f) for f in chi]
    dlep=np.zeros_like(lep)
    cm_nodes,cm_weights=leggauss(nc)
    cosphi=np.cos(2*math.pi*np.arange(nphi)/nphi)[None,:]
    mu_nodes,mu_weights=leggauss(nmu)
    entropy_direct=0.
    for species,(M,grid,ff) in enumerate(zip(cfg["masses"],cfg["chi_grid"],chi)):
      fch=interp_fermi(grid,ff)
      px,wp=radial_nodes(M,quadrature_cut,nr,scheme=radial_scheme,segments=segments)
      kx,wk=radial_nodes(cfg["mlepton"],quadrature_cut,nr,scheme=radial_scheme,segments=segments)
      for p,ap in zip(px,wp):
        e=math.hypot(p,M)
        for k,ak in zip(kx,wk):
          o=math.hypot(k,cfg["mlepton"])
          for mu,au in zip(mu_nodes,mu_weights):
            tot=e+o
            P=math.sqrt(p*p+k*k+2*p*k*mu)
            s=tot*tot-P*P
            ps2=(s-(M+cfg["mlepton"])**2)*(s-(M-cfg["mlepton"])**2)/(4*s)
            if ps2<=0:continue
            ps=math.sqrt(ps2)
            beta=P/tot
            gamma=tot/math.sqrt(s)
            proj=p*(p+k*mu)/P
            ca=max(-1.,min(1.,gamma*(proj-beta*e)/ps))
            estar=(s+M*M-cfg["mlepton"]**2)/(2*math.sqrt(s))
            A=2*ps2/mv**2
            ymax=math.log1p(2*A)
            y=(cm_nodes+1)*ymax/2
            c=(1-np.expm1(y)/A)[:,None]
            wcm=(cm_weights*ymax/(2*A)*np.exp(-y))[:,None]
            ep=gamma*(estar+beta*ps*(ca*c+math.sqrt(max(0,1-ca*ca))*
                                            np.sqrt(np.maximum(0,1-c*c))*cosphi))
            op=tot-ep
            t=-2*ps2*(1-c)
            D=s-M*M-cfg["mlepton"]**2
            U=-D-t
            amp_num=2*G*G*(D*D+U*U+2*(M*M+cfg["mlepton"]**2)*t)
            fi=float(fch(e))
            lj=float(lepf(o))
            fo=fch(ep)
            lo=lepf(op)
            F=fi*lj*(1-fo)*(1-lo)
            R=fo*lo*(1-fi)*(1-lj)
            J=F-R
            v=math.sqrt(max(0.,(e*o-p*k*mu)**2-M*M*cfg["mlepton"]**2))/(e*o)
            factor=16/(512*math.pi**5)*ap*p*p*ak*k*k*au*v/s
            event_weight=factor*wcm*amp_num/nphi
            rates=np.broadcast_to(event_weight,ep.shape)*J
            deposit(dchi[species],grid,ep,rates.ravel())
            deposit(dchi[species],grid,[e],[-float(np.sum(rates))])
            deposit(dlep,cfg["lg"],op,rates.ravel())
            deposit(dlep,cfg["lg"],[o],[-float(np.sum(rates))])
            entropy_direct+=float(np.sum(rates*np.log(F/R)))
    dfdchi=[x/w for x,w in zip(dchi,cfg["chi_w"])]
    dfdlep=dlep/cfg["lw"]
    entropy_node=-sum(float(np.dot(dn,logit(f))) for dn,f in zip(dchi,chi))
    entropy_node-=float(np.dot(dlep,logit(lep)))
    return dfdchi,dfdlep,{"entropy_direct":entropy_direct,
        "entropy_node":entropy_node,
        "Qchi":[float(np.dot(g,dn)) for g,dn in zip(cfg["chi_grid"],dchi)],
        "Qlep":float(np.dot(cfg["lg"],dlep)),
        "Nchi_rates":[float(np.sum(dn)) for dn in dchi],
        "Nlep_rate":float(np.sum(dlep))}


def advance(cfg,chi,lep,dt):
    dc,dl,_=rhs(cfg,chi,lep)
    cnext=[f+dt*d for f,d in zip(chi,dc)]
    lnext=lep+dt*dl
    if any(np.any((f<=0)|(f>=1)) for f in [*cnext,lnext]):
        return None
    if ledger(cfg,cnext,lnext)["entropy"]<ledger(cfg,chi,lep)["entropy"]-1e-12:
        return None
    return cnext,lnext


def evolve(cfg,total_time=1e10,atol=1e-7,dt_initial=1e10,max_steps=250):
    chi=[x.copy() for x in cfg["initial_chi"]]
    lep=cfg["initial_lepton"].copy()
    start=ledger(cfg,chi,lep)
    t=0.;dt=dt_initial;accepted=rejected=0
    smallest_entropy=float("inf");spectral_error=0.
    for _ in range(max_steps):
      if t>=total_time*(1-1e-14):break
      dt=min(dt,total_time-t)
      full=advance(cfg,chi,lep,dt)
      half=advance(cfg,chi,lep,dt/2)
      pair=advance(cfg,*half,dt/2) if half is not None else None
      if full is not None and pair is not None:
        err=max(max(float(np.max(abs(f-g))) for f,g in zip(full[0],pair[0])),
                float(np.max(abs(full[1]-pair[1]))))
        if err<=atol:
          sb=ledger(cfg,chi,lep)["entropy"]
          chi,lep=pair
          sa=ledger(cfg,chi,lep)["entropy"]
          smallest_entropy=min(smallest_entropy,sa-sb)
          spectral_error=max(spectral_error,err)
          t+=dt;accepted+=1
          dt*=min(2.,max(1.1,.9*math.sqrt(atol/max(err,1e-30))))
          continue
      rejected+=1;dt*=.5
      if dt<1e-18:raise ArithmeticError("positive entropy step unavailable")
    else:raise RuntimeError("finite adaptive steps exhausted")
    end=ledger(cfg,chi,lep)
    return {"accepted":accepted,"rejected":rejected,"time":t,
        "min_step_entropy_gain":smallest_entropy,
        "max_local_spectral_error":spectral_error,
        "Nchi1_drift":end["Nchi1"]-start["Nchi1"],
        "Nchi2_drift":end["Nchi2"]-start["Nchi2"],
        "Nlepton_drift":end["Nlepton"]-start["Nlepton"],
        "energy_drift":end["Etotal"]-start["Etotal"],
        "entropy_gain":end["entropy"]-start["entropy"],
        "final_fchi1_min":float(min(chi[0])),
        "final_fchi2_min":float(min(chi[1])),
        "final_lepton_min":float(min(lep)),
        "not_real_phi_or_frw":True}
