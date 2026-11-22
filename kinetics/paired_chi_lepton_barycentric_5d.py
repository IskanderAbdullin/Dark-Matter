"""Conservative χ+lepton 5D WEAK-to-binned-PROJECTED event rates.

Use exact same on-shell finite-vector 2->2 microevents and Pauli flux
as nonthermal_chi_lepton_5d_weak. Deposit source/sink of χ and lepton
onto separate physical energy grids by nonnegative barycentric weights.
Partition unity and first energy moment EXACT within the truncated
incoming integral, so weak Qχ=-Qℓ matches projected strong-rate moments.

CRITICAL LIMIT: evaluated Fermi occupations remain prescribed analytic
f(E), not the bin populations. The projected dN_i are NOT a validated
nonlinear mode-resolved strong Cχ/Cℓ; barycentric bin-interpolated
entropy can VIOLATE the H theorem for nonthermal curved log odds,
despite exact number and energy. Do not advance these rates in time
as physical Pauli kinetics without an entropy-consistent projection.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.nonthermal_chi_lepton_5d_weak import event,gauss


def deposit_linear(rate_bins,energy_grid,energies,amounts):
    """Positive partition unity and exact first-energy-moment weights."""
    energy=np.atleast_1d(np.asarray(energies,float))
    rate=np.atleast_1d(np.asarray(amounts,float))
    if energy.shape!=rate.shape or not np.all(np.isfinite(energy+rate)):
        raise ValueError("matching finite energy and event rate vectors")
    lo=float(energy_grid[0]);hi=float(energy_grid[-1])
    width=float(energy_grid[1]-lo)
    if np.any(energy<lo-1e-10) or np.any(energy>hi+1e-10):
        raise ValueError("physical incoming/outgoing energy outside projection domain")
    z=(energy-lo)/width
    z=np.maximum(0.,np.minimum(z,len(energy_grid)-1-1e-13))
    index=np.floor(z).astype(int)
    frac=z-index
    np.add.at(rate_bins,index,(1-frac)*rate)
    np.add.at(rate_bins,index+1,frac*rate)


def projected_pair(fchi,flepton,*,M=.1,m=0.,Te=1.,G=.001,mV=.1,
                   nr=4,nmu=4,nphi=12,nang=12,cut_kinetic_T=12.,
                   bins=20):
    """Integrated events and conservative χ/lepton energy-binned source.

    NOTE: scalar bins are *total weighted particle-number rates*, not
    physical occupation df/dt without density-of-states and nonlinear
    Fermi-aware deposition. Width up to total available incoming energy.
    """
    if not(M>0 and m>=0 and Te>0 and G>0 and mV>0 and nr>=4
           and nmu>=4 and nphi>=8 and nang>=4 and cut_kinetic_T>0
           and isinstance(bins,int) and bins>=4):
        raise ValueError("physical inputs and bin count")
    max_total=M+m+2*cut_kinetic_T*Te
    gridchi=np.linspace(M,max_total,bins+1)
    gridlep=np.linspace(m,max_total,bins+1)
    nchi=np.zeros(bins+1)
    nlep=np.zeros(bins+1)
    pmax=math.sqrt((M+cut_kinetic_T*Te)**2-M*M)
    kmax=math.sqrt((m+cut_kinetic_T*Te)**2-m*m)
    ps,wp=gauss(0.,pmax,nr)
    ks,wk=gauss(0.,kmax,nr)
    mus,wmu=gauss(-1.,1.,nmu)
    cos,wcos=leggauss(nang)
    phi=2*math.pi*np.arange(nphi)/nphi
    pref=16/(512*math.pi**5)
    weak_heat=weak_entropy=0.
    for p,ap in zip(ps,wp):
        for k,ak in zip(ks,wk):
            for u,au in zip(mus,wmu):
                for c,ac in zip(cos,wcos):
                    ev=event(M,m,float(p),float(k),float(u),G,mV,
                             float(c),phi,fchi,flepton)
                    weight=pref*ap*p*p*ak*k*k*au*ac*ev["vm"]/ev["s"]
                    event_rate=weight*ev["amplitude"]*ev["flux"]/nphi
                    deposit_linear(nchi,gridchi,ev["Ep"],event_rate)
                    deposit_linear(nchi,gridchi,[ev["E"]],
                                   [-float(np.sum(event_rate))])
                    deposit_linear(nlep,gridlep,ev["omega_p"],event_rate)
                    deposit_linear(nlep,gridlep,[ev["omega"]],
                                   [-float(np.sum(event_rate))])
                    weak_heat+=float(np.dot(event_rate,ev["dEchi"]))
                    weak_entropy+=weight*ev["amplitude"]*float(
                        np.mean(ev["entropy"]))
    fchi_bins=np.asarray(fchi(gridchi),float)
    flep_bins=np.asarray(flepton(gridlep),float)
    if (np.any(fchi_bins<=0) or np.any(fchi_bins>=1) or
        np.any(flep_bins<=0) or np.any(flep_bins>=1)):
        raise ValueError("strict interior Fermi statistics required for entropy")
    strong_entropy=(
        -float(np.dot(nchi,np.log(fchi_bins/(1-fchi_bins))))
        -float(np.dot(nlep,np.log(flep_bins/(1-flep_bins)))))
    return {
      "Echi_grid":gridchi,"Elepton_grid":gridlep,
      "dNchi_bins":nchi,"dNlepton_bins":nlep,
      "number_chi_residual":float(np.sum(nchi)),
      "number_lepton_residual":float(np.sum(nlep)),
      "weak_Qchi":float(weak_heat),
      "projected_Qchi":float(np.dot(nchi,gridchi)),
      "projected_Qlepton":float(np.dot(nlep,gridlep)),
      "direct_event_entropy":float(weak_entropy),
      "interpolated_bin_entropy":float(strong_entropy),
      "entropy_projection_relative_error":
          float(strong_entropy/weak_entropy-1)
          if weak_entropy!=0 else None,
      "not_a_nonlinear_Pauli_occupation_strong_operator":True,
      "negative_H_theorem_possible_on_curved_nonthermal_log_odds":True,
      "class_ready":False
    }


if __name__=="__main__":
    import json
    from scipy.special import expit
    def pack(x):
        return {k:v for k,v in x.items() if not isinstance(v,np.ndarray)}
    smoothchi=lambda E:expit(-np.asarray(E)/.9+.15*np.exp(-np.asarray(E)/.9))
    smoothlep=lambda E:expit(-np.asarray(E)+.05*np.exp(-np.asarray(E)))
    oscillatorychi=lambda E:expit(-np.asarray(E)+.3*np.sin(4*np.asarray(E)))
    oscillatorylep=lambda E:expit(-np.asarray(E)+.3*np.cos(4*np.asarray(E)))
    print(json.dumps({
       "smooth":{str(b):pack(projected_pair(smoothchi,smoothlep,bins=b))
                 for b in (10,20,40,80,160)},
       "oscillatory_negative_control":{
          str(b):pack(projected_pair(oscillatorychi,oscillatorylep,bins=b))
          for b in (10,20,40,80,160,320)},
       "not_full_nonthermal_chi_lepton_strong_collision":True
    },indent=2))
