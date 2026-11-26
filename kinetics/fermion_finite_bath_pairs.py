"""Conservative implicit Pauli pair kinetics with finite radiation reservoir.

REDUCED two-state/network benchmark ONLY, NOT microscopic 5D scattering.
For fermion energy levels E_i<E_j with number capacities V_i,V_j:
J=Gamma exp[-(E_j-E_i)/T] f_i(1-f_j)-Gamma f_j(1-f_i).
Transferred number dN solves dN=dt J at NEW fermion occupations
and NEW finite bath T=((E_SM-dN*DeltaE)/A)^(1/4). The scalar Brent
bracket enforces Pauli 0<=f<=1 and E_SM>=0 at arbitrarily large dt.
Total fermion number and dark+SM energy are conserved to roundoff.
For radiation E_SM=A*T^4, S_SM=4/3 A*T^3 at fixed scale factor, entropy
production J ln(r_up/r_down)>=0 under the pair dynamics.

Gamma here is an arbitrary positive pair coefficient, NOT CALIBRATED
from actual χe 5D invariant collision and NOT a real early-Universe SM
or Bose dark-matter production model. Reduced network is a numerical
conservation and integration reference, not a physical gamma_el.
"""
from __future__ import annotations
import math
from scipy.optimize import brentq


def _fermi_entropy(f):
    if f<=0 or f>=1:
        return 0.
    return -f*math.log(f)-(1-f)*math.log1p(-f)


def entropy_total(f,capacities,ESM,A):
    if ESM<0:
        return float("nan")
    return (sum(V*_fermi_entropy(x) for x,V in zip(f,capacities))
            +4*A/3*(ESM/A)**.75)


def evolve_pair(Ei,Ej,fi,fj,Vi,Vj,ESM,A,Gamma,dt):
    vals=(Ei,Ej,fi,fj,Vi,Vj,ESM,A,Gamma,dt)
    if (not all(math.isfinite(x) for x in vals) or not 0<=Ei<Ej or
            not 0<=fi<=1 or not 0<=fj<=1 or
            Vi<=0 or Vj<=0 or ESM<=0 or A<=0 or Gamma<0 or dt<0):
        raise ValueError("Finite ordered energies, Pauli bounds, positive bath/capacity, nonnegative rate/dt")
    deltaE=Ej-Ei
    E_dark0=Vi*fi*Ei+Vj*fj*Ej
    S0=entropy_total((fi,fj),(Vi,Vj),ESM,A)
    if Gamma==0 or dt==0:
        return dict(fi=fi,fj=fj,ESM=ESM,T=(ESM/A)**.25,dN=0.,
         particle_residual=0.,energy_residual=0.,entropy_delta=0.,
         balance_residual=0.,reference="reduced network, not microscopic 5D collision")
    lo=-min(Vj*fj,Vi*(1-fi))
    hi=min(Vi*fi,Vj*(1-fj),ESM/deltaE)
    def states(dN):
        ni=fi-dN/Vi
        nj=fj+dN/Vj
        reservoir=ESM-dN*deltaE
        T=(max(0.,reservoir)/A)**.25
        upward=math.exp(-deltaE/T) if T>0 else 0.
        flow=Gamma*(upward*ni*(1-nj)-nj*(1-ni))
        return ni,nj,reservoir,T,flow
    def residual(dN):
        return dN-dt*states(dN)[-1]
    transfer=brentq(residual,lo,hi,xtol=5e-15,rtol=5e-15)
    ni,nj,newE,newT,flow=states(transfer)
    if (not -1e-13<=ni<=1+1e-13 or not -1e-13<=nj<=1+1e-13
            or newE<-1e-12):
        raise ArithmeticError("Unphysical implicit pair solution")
    E_dark1=Vi*ni*Ei+Vj*nj*Ej
    Nerr=Vi*ni+Vj*nj-Vi*fi-Vj*fj
    Eerr=E_dark1+newE-E_dark0-ESM
    Snew=entropy_total((ni,nj),(Vi,Vj),max(0.,newE),A)
    return dict(fi=ni,fj=nj,ESM=newE,T=newT,dN=transfer,
      particle_residual=Nerr,energy_residual=Eerr,
      entropy_delta=Snew-S0,balance_residual=transfer-dt*flow,
      reference="reduced finite-bath Pauli pair, not full 5D microphysics")


def evolve_network(E,capacities,f,ESM,A,edges,dt,steps=1):
    """Lie-split implicit rate network; prescribed transitions are not physical collision kernel.

    Each edge=(i,j,Gamma), i<j. Order effects must be tested
    independently from energy/number/entropy preservation.
    """
    E=list(map(float,E));V=list(map(float,capacities))
    pop=list(map(float,f))
    if (len(E)<2 or len(E)!=len(V) or len(E)!=len(pop) or
            any(E[i+1]<=E[i] for i in range(len(E)-1)) or
            any(v<=0 for v in V) or any(not 0<=v<=1 for v in pop) or
            type(steps) is not int or steps<1 or dt<0):
        raise ValueError("Increasing energy grid, matching Pauli capacities and positive step")
    tot0=sum(v*x for v,x in zip(V,pop))
    en0=sum(v*x*e for v,x,e in zip(V,pop,E))+ESM
    sn0=entropy_total(pop,V,ESM,A)
    for _ in range(steps):
        for i,j,gamma in edges:
            if not 0<=i<j<len(E):
                raise ValueError("ordered edge i<j required")
            r=evolve_pair(E[i],E[j],pop[i],pop[j],V[i],V[j],
                          ESM,A,gamma,dt/steps)
            pop[i]=r["fi"];pop[j]=r["fj"];ESM=r["ESM"]
    tot=sum(v*x for v,x in zip(V,pop))
    energy=sum(v*x*e for v,x,e in zip(V,pop,E))+ESM
    return dict(f=pop,ESM=ESM,T=(ESM/A)**.25,
        number_residual=tot-tot0,total_energy_residual=energy-en0,
        entropy_delta=entropy_total(pop,V,ESM,A)-sn0,
        reference="reduced network with arbitrary pair rates, not microscopic collision integral")
