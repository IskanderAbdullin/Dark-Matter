"""Conservative weak-event energy-bin projection; exact FD Fermi logit interpolation.

This module is NOT a physical finite-volume time integrator.
For each event rate W, deposit W*[B_j(Eout)-B_j(Ein)] using a linear
partition of unity in physical E.  The discrete number and energy
moments match event-number and physical event-energy exactly, with
SM heat set to minus dark energy; explicit Euler may still violate
Pauli and a projected source is not automatically an H-theorem.

Separately, interpolate log(f/(1-f)) linearly in physical energy.
An exact FD has logit(f)=(mu-E)/T, and hence remains EXACT under any
properly sampled grid.  Unlike PCHIP(f(q)), this does not create a
spurious thermal collision source from pure FD curvature on coarse q.
Requires strictly 0<f<1 and explicit tail-extrapolation assumptions.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit,logit


def logit_energy_interpolator(q,f,mass):
    """Pauli-bounded interpolator, exact for arbitrary FD(T,mu) in E.

    Strictly excludes endpoints f=0,1 until their limiting treatment
    is derived. Extrapolating the endpoint logit slopes is an assumption
    for a NONTHERMAL high-momentum distribution, not a physical result.
    """
    q=np.asarray(q,dtype=float)
    f=np.asarray(f,dtype=float)
    if (q.ndim!=1 or len(q)<3 or f.shape!=q.shape or
        not np.isfinite(mass) or mass<=0 or
        np.any(~np.isfinite(q)) or np.any(~np.isfinite(f)) or
        q[0]!=0 or np.any(np.diff(q)<=0) or
        np.any((f<=0)|(f>=1))):
        raise ValueError("strict finite 0<f<1, q0=0, increasing grid, M>0 required")
    energy=np.hypot(q,mass)
    odds=logit(f)
    def evaluate(momentum):
        p=np.asarray(momentum,float)
        if np.any(~np.isfinite(p)) or np.any(p<0):
            raise ValueError("finite nonnegative momentum required")
        e=np.hypot(p,mass)
        ix=np.searchsorted(energy,e,side="right")-1
        ix=np.clip(ix,0,len(energy)-2)
        interpolated=odds[ix]+(odds[ix+1]-odds[ix])*(
            e-energy[ix])/(energy[ix+1]-energy[ix])
        return expit(interpolated)
    return evaluate


def _energy_weights(E,nodes):
    if not np.isfinite(E) or E<nodes[0] or E>nodes[-1]:
        raise ValueError("event energy outside conservative grid; extend grid")
    i=np.searchsorted(nodes,E,side="right")-1
    i=int(np.clip(i,0,len(nodes)-2))
    frac=float((E-nodes[i])/(nodes[i+1]-nodes[i]))
    w=np.zeros_like(nodes)
    w[i]=1-frac
    w[i+1]=frac
    return w


def project_weak_events_energy_bins(q_nodes,mass,energies_initial,
                                    energies_final,number_rates):
    """A conservative event *moment projection*, NOT stable time evolution.

    The passed signed event number rates already include physical
    amplitudes/quantum factors if used for real scattering.  This module
    does not compute them. It fails closed if any event is outside grid.
    """
    q=np.asarray(q_nodes,float)
    if (q.ndim!=1 or len(q)<3 or np.any(~np.isfinite(q)) or
        q[0]!=0 or np.any(np.diff(q)<=0) or
        not np.isfinite(mass) or mass<=0):
        raise ValueError("valid q nodes starting zero, mass>0 required")
    egrid=np.hypot(q,mass)
    Ei=np.asarray(energies_initial,float)
    Ef=np.asarray(energies_final,float)
    W=np.asarray(number_rates,float)
    if (Ei.ndim!=1 or Ei.shape!=Ef.shape or Ei.shape!=W.shape or
        np.any(~np.isfinite(Ei)) or np.any(~np.isfinite(Ef)) or
        np.any(~np.isfinite(W))):
        raise ValueError("finite matching 1D energy/rate arrays required")
    dN=np.zeros_like(q)
    for Ein,Eout,rate in zip(Ei,Ef,W):
        dN+=rate*(_energy_weights(float(Eout),egrid)-
                 _energy_weights(float(Ein),egrid))
    expected=float(np.dot(W,Ef-Ei))
    measured=float(np.dot(egrid,dN))
    return {
      "source_cell_number_per_time":dN,
      "number_residual":float(np.sum(dN)),
      "dark_energy_rate_from_grid":measured,
      "dark_energy_rate_from_exact_events":expected,
      "sm_energy_rate":-measured,
      "energy_projection_defect":measured-expected,
      "not_time_integrator":True,
      "not_Pauli_stable_time_update":True,
    }
