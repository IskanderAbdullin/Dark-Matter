"""Physical full-Pauli 5D χe equilibrium event conductances projected to energy bins.

NOT an exact nonlinear finite-bin Boltzmann operator. Exact physical recoil
E'_chi-E_chi is deposited by barycentric weights; that preserves mean energy
but ADDS positive Jensen variance to the energy-mode Dirichlet response.
Raw projected equilibrium transitions are NOT microscopically reversible
on finite bins. A symmetric closure enforces discrete FD detailed balance,
but alters the physical kernel; both defects are measured, never hidden.

All energies GeV, G GeV^-2, conductances and D GeV^4. Broken-EW vector
chi-e massless electron contact only, no actual SM/FRW or 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit
from numpy.polynomial.legendre import leggauss


def gauss(n, a, b):
    x,w=leggauss(n)
    return (a+b)/2+(b-a)*x/2, (b-a)*w/2


def cm_event(M,p,k,mu,cost,phi):
    E=math.hypot(M,p)
    rx=k*math.sqrt(max(0.,1-mu*mu)); rz=p+k*mu; Et=E+k
    s=M*M+2*k*(E-p*mu); root=math.sqrt(s)
    ps=(s-M*M)/(2*root); Es=(s+M*M)/(2*root)
    factor=p*rz/(root*(Et+root))-E/root
    cx=factor*rx; cz=p+factor*rz; clen=math.hypot(cx,cz)
    if clen<=0: raise ArithmeticError("degenerate CM orientation")
    nx,nz=cx/clen,cz/clen
    ct=np.asarray(cost); ang=np.asarray(phi)
    st=np.sqrt(np.maximum(0,1-ct*ct))
    xp=ps*(ct*nx-st*np.cos(ang)*nz)
    zp=ps*(ct*nz+st*np.cos(ang)*nx)
    dot=rx*xp+rz*zp
    Ep=(Et*Es+dot)/root
    return Ep,Et-Ep,-2*ps*ps*(1-ct)


def graph(*,M=.1,T=1.,G=.001,gchi=4.,ge=4.,nchi=32,pmax=22.,
          nk=12,nmu=10,ncm=10,nphi=10,kmax=22.,stretch=1.6):
    """Integrate true 5D vector-contact |M|² over incident p cells and e FD.

    The raw directed conductance matrix is a POSITIVE coarse quadrature
    of F_eq times exact scattering, before any reversible closure.
    """
    scalars=(M,T,G,gchi,ge,pmax,kmax,stretch)
    if (not all(math.isfinite(x) for x in scalars) or M<=0 or T<=0 or G<0
        or gchi<=0 or ge<=0 or pmax<=0 or kmax<=0 or stretch<=0
        or any(type(n)!=int or n<3 for n in (nchi,nk,nmu,ncm,nphi))):
        raise ValueError("finite physical inputs and grids >=3 required")
    p=pmax*(np.arange(nchi)/(nchi-1))**stretch
    E=np.hypot(M,p)
    edges=np.concatenate(([0.],(p[:-1]+p[1:])/2,[pmax]))
    # Exact state capacities for piecewise radial p^3 bins.
    caps=gchi/(6*math.pi**2)*(edges[1:]**3-edges[:-1]**3)
    ks,kw=gauss(nk,0,kmax)
    mus,mw=gauss(nmu,-1,1)
    us,uw=gauss(ncm,-1,1)
    phis=2*np.pi*(np.arange(nphi)+.5)/nphi
    U=us[:,None];Phi=phis[None,:]
    angle_weights=uw[:,None]*np.ones((1,nphi))/nphi
    L=np.zeros((nchi,nchi),float)
    physical_D_all=physical_D_inside=outweight=allweight=0.
    for i,(pi,Ei,Vi) in enumerate(zip(p,E,caps)):
        fi=float(expit(-Ei/T))
        for k,wk in zip(ks,kw):
            fe=float(expit(-k/T))
            for mu,wmu in zip(mus,mw):
                s=M*M+2*k*(Ei-pi*mu)
                flux=1-pi*mu/Ei
                Ep,kp,t=cm_event(M,float(pi),float(k),float(mu),U,Phi)
                Delta=s-M*M; usmand=2*M*M-s-t
                amp=2*G*G*(Delta*Delta+(usmand-M*M)**2+2*M*M*t)
                dsdom=amp/(64*math.pi**2*s)
                Feq=fi*fe*(1-expit(-Ep/T))*(1-expit(-kp/T))
                # V_i includes gchi; ge/(2pi) matches previous strong Cchi.
                W=Vi*ge/(2*math.pi)*wk*k*k*wmu*flux*dsdom*Feq*angle_weights
                evE=np.broadcast_to(Ep,W.shape).ravel()
                evW=W.ravel()
                physical_D_all+=.5*float(np.dot(evW,((evE-Ei)/T)**2))
                allweight+=float(np.sum(evW))
                inside=(evE<=E[-1])&(evE>=E[0]-1e-10)
                outweight+=float(np.sum(evW[~inside]))
                evE=np.clip(evE[inside],E[0],E[-1]);evW=evW[inside]
                physical_D_inside+=.5*float(np.dot(evW,((evE-Ei)/T)**2))
                j=np.searchsorted(E,evE,side="right")-1
                j=np.clip(j,0,nchi-2)
                a=(evE-E[j])/(E[j+1]-E[j])
                np.add.at(L,(i,j),evW*(1-a))
                np.add.at(L,(i,j+1),evW*a)
    sym=(L+L.T)/2
    d=(E[None,:]-E[:,None])/T
    Dmap=.5*float(np.sum(L*d*d))
    total=float(np.sum(L))
    net=np.sum(L.T-L,axis=1)
    raw_heat=float(np.sum(L*(E[None,:]-E[:,None])))
    if allweight==0:
        return {"E":E,"p":p,"capacities":caps,"directed":L,"symmetric":sym,
                "D_true_all":0.,"D_true_inside":0.,"D_graph":0.,
                "D_graph_symmetric":0.,"jensen_excess":0.,
                "jensen_relative":0.,"unpaired_asymmetry":0.,
                "raw_stationary_abs_number_source_GeV4":0.,
                "raw_stationary_relative_L1":0.,
                "raw_FD_heat_without_reversible_closure_GeV5":0.,
                "outgoing_out_of_support_weight_fraction":0.,
                "accepted_forward_weight_GeV4":0.,"all_forward_weight_GeV4":0.,
                "params":dict(M=M,T=T,G=G,nchi=nchi,pmax=pmax,nk=nk,
                              nmu=nmu,ncm=ncm,nphi=nphi,kmax=kmax,stretch=stretch),
                "is_exact_nonlinear_microphysical_bin_kernel":False}
    return {
        "p":p,"E":E,"capacities":caps,"directed":L,"symmetric":sym,
        "D_true_all":physical_D_all,"D_true_inside":physical_D_inside,
        "D_graph":Dmap,"D_graph_symmetric":.5*float(np.sum(sym*d*d)),
        "jensen_excess":Dmap-physical_D_inside,
        "jensen_relative":(Dmap-physical_D_inside)/physical_D_inside,
        "unpaired_asymmetry":float(np.sum(np.abs(L-L.T)))/(2*total),
        "raw_stationary_abs_number_source_GeV4":float(np.sum(np.abs(net))),
        "raw_stationary_relative_L1":float(np.sum(np.abs(net)))/total,
        "raw_FD_heat_without_reversible_closure_GeV5":raw_heat,
        "outgoing_out_of_support_weight_fraction":outweight/allweight,
        "accepted_forward_weight_GeV4":total,
        "all_forward_weight_GeV4":allweight,
        "params":dict(M=M,T=T,G=G,nchi=nchi,pmax=pmax,nk=nk,
                      nmu=nmu,ncm=ncm,nphi=nphi,kmax=kmax,stretch=stretch),
        "is_exact_nonlinear_microphysical_bin_kernel":False,
    }


def calibrated_pair_rates(g):
    """Rates derived from symmetrized equilibrium L, at one T only.

    Such Γ are NOT directly computed nonlinear 5D χe rates at other T/fχ.
    """
    E=g["E"];L=g["symmetric"];T=g["params"]["T"]
    f=expit(-E/T)
    return [(i,j,float(L[i,j]/(f[j]*(1-f[i]))))
            for i in range(len(E)) for j in range(i+1,len(E)) if L[i,j]>0]


def graph_instantaneous_heat(g,f):
    """Reversibly closed pair approximation; NOT exact nonthermal Cchi."""
    E=g["E"];T=g["params"]["T"]
    f=np.asarray(f,float)
    if f.shape!=E.shape or np.any(~np.isfinite(f)) or np.any((f<0)|(f>1)):
        raise ValueError("invalid Fermi occupation")
    Q=0.;net=np.zeros_like(f)
    for i,j,Gamma in calibrated_pair_rates(g):
        up=Gamma*math.exp(-(E[j]-E[i])/T)*f[i]*(1-f[j])
        down=Gamma*f[j]*(1-f[i])
        J=up-down
        net[i]-=J;net[j]+=J
        Q+=(E[j]-E[i])*J
    return Q,net


def local_pair_implicit_step(g,i,j,fi,fj,esm,radiation_A,dt):
    """Finite bath with microscopic equilibrium Γ FROZEN at reference T.

    Uses independently tested implicit pair reference; bath response outside
    reference T is not rederived from the changing electron distribution.
    """
    from fermion_finite_bath_pairs import evolve_pair
    if not 0<=i<j<len(g["E"]): raise ValueError("ordered pair required")
    E=g["E"];f=expit(-E/g["params"]["T"]);L=g["symmetric"][i,j]
    Gamma=float(L/(f[j]*(1-f[i]))) if L>0 else 0.
    result=evolve_pair(float(E[i]),float(E[j]),fi,fj,
                       float(g["capacities"][i]),float(g["capacities"][j]),
                       esm,radiation_A,Gamma,dt)
    result["reference_temperature_GeV"]=g["params"]["T"]
    result["frozen_microphysical_linear_conductance_GeV4"]=float(L)
    result["Gamma_at_reference_T"]=Gamma
    result["physical_nonthermal_full_collision"]=False
    return result
