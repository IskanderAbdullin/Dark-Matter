"""Test whether energy and a light cap alone make a high inverse source small."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from three_species_decay_onsager_reference import M, P, m, mp
from warm_light_occupation_stage_audit_reference import read


def radial(center,width,mass):
    # Retain the centered width explicitly even if a float would round an
    # endpoint to the center at very large light momentum.
    x,w=roots_legendre(16)
    E=center+.5*width*x
    return .5*width*np.sum(w*E**2*np.sqrt((E-mass)*(E+mass)))/(2*np.pi**2)


def one(Ea,Eparent,high):
    D=M*M-m*m-mp*mp
    center_s=D/Eparent/2
    center_l=Eparent-center_s
    ws=.12*center_s
    # Commit only 1/4 of the total energy allowance to each fixed bath.
    incoming=Ea/4
    delta=.5
    def f(logw):
        return np.log(4*delta*radial(center_l,np.exp(logw),m)/incoming)
    guess=incoming/((4*delta)/(2*np.pi**2)*center_l**3)
    wl=float(np.exp(brentq(f,np.log(guess)-4,np.log(guess)+4,xtol=4e-14)))
    fs=incoming/radial(center_s,ws,mp)
    assert abs(4*delta*radial(center_l,wl,m)/incoming-1.)<5e-13
    lo_l,hi_l=np.nextafter(center_l-.5*wl,-np.inf),np.nextafter(center_l+.5*wl,np.inf)
    lo_s,hi_s=np.nextafter(center_s-.5*ws,-np.inf),np.nextafter(center_s+.5*ws,np.inf)
    Elo,Ehi=lo_l+lo_s,hi_l+hi_s
    assert Elo>np.hypot(M,high)
    # Check a conservative admissible rectangle for the two daughter
    # energies. The max of the invariant minimum is not assumed at a
    # corner: use the upper energy product and lower momentum product.
    kl=np.sqrt((lo_l-m)*(lo_l+m))
    ks=np.sqrt((lo_s-mp)*(lo_s+mp))
    margin=(D/2-(hi_l*hi_s-kl*ks),lo_l*lo_s+kl*ks-D/2)
    assert min(margin)>0.
    # In u=E_light and e=E_scalar, the positive inverse source has
    # measure M/(pi² P) * (u+e) after division by the rest width times M.
    Q=M/(np.pi**2*P)*delta*fs*wl*ws*Eparent
    x,w=roots_legendre(16)
    sum_energy=center_l+center_s+.5*(wl*x[:,None]+ws*x[None,:])
    cross_check=float(M/(np.pi**2*P)*delta*fs*.25*wl*ws*(w@sum_energy@w))
    assert abs(cross_check/Q-1.)<6e-15
    # Gamma M / E times (1 + <b>) is a loss bound for the frozen
    # combined affine operator. Use whole scalar width even at those
    # parent momenta which would miss a fraction of the scalar band.
    ph=np.sqrt((Elo-M)*(Elo+M))
    B=M/(2*ph*P)*fs*ws
    coeff=(1+B)/Elo
    lo_stationary=Q/coeff/Ea
    return dict(parent_energy_center_GeV=Eparent,parent_band_GeV=[Elo,Ehi],
                light_and_scalar_energy_centers_GeV=[center_l,center_s],
                light_and_scalar_widths_GeV=[wl,ws],
                relative_allocations=[.25,.25],scalar_per_mode=fs,light_per_mode=delta,
                invariant_half_gap_GeV2=margin,positive_energy_input_over_GammaM=Q,
                scalar_angle_average_ceiling=B,max_decay_coefficient_over_GammaM=coeff,
                frozen_reservoir_lower_over_Ea=lo_stationary)


def run(source):
    a=read(source)["branches"][-1]
    assert a["label"]=="q64-24"
    E=a["energy_upper_GeV4"]
    cut=1./a["entropy_expansion_from_original_anchor"]
    values=[one(E,h,cut) for h in (.75,1.,2.,10.,30.,100.,200.,300.)]
    return dict(scope="Synthetic fixed daughter reservoirs consuming one half of the late gross budget. "
                      "The implied stationary heavy reservoir ignores backreaction and need not obey the "
                      "same total budget. This cannot be used as a physical heavy abundance or allowed tail. "
                      "It tests only whether a light cap and energy moments alone enforce a small inverse source.",
                snapshot_sha256=a["snapshot_sha256"],energy_allowance_GeV4=E,
                heavy_momentum_cut_at_late_anchor_GeV=cut,controls=values)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.source),indent=2,allow_nan=False))
