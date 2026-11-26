"""Conditional neutrino inverse-decay vs H diagnostic, NOT actual quantum rate.

For a broken-phase Z' with two active ν flavors g', vacuum rest
Gamma_nunu=gprime^2*mV/(12π). The *Maxwell-Boltzmann* thermal
time-dilated decay is Gamma_rest*K1(mV/T)/K2(mV/T);
this omits final ν Pauli, inverse νν→V and vector Bose and cannot
prove equilibration or cosmological Neff. Radiation H uses assumed gstar.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import kv

ME_GEV=.00051099895
MPL_GEV=1.2209e19


def boltzmann_time_dilation(mV,T,*,n=256):
    """Independent p/T integral of <mV/E> over MB massive vector."""
    if mV<=0 or T<=0: raise ValueError
    z,w=leggauss(n);q=30*(z+1);w=30*w
    x=mV/T
    eps=np.hypot(q,x)
    ww=w*q*q*np.exp(-eps)
    return float(np.dot(ww,x/eps)/np.sum(ww))


def analytic_time_dilation(mV,T):
    x=mV/T
    if not (x>0):raise ValueError
    return float(kv(1,x)/kv(2,x))


def conditional_neutrino_decay(*,mV_me=.03,T_me=1.,
                                G_me_minus2=.001,gchi_eff=None,
                                gstar=10.75):
    if not all(math.isfinite(x) and x>0 for x in
               (mV_me,T_me,G_me_minus2,gstar)):
        raise ValueError
    if gchi_eff is None:gchi_eff=math.sqrt(4*math.pi)
    if not math.isfinite(gchi_eff) or gchi_eff<=0:
        raise ValueError
    mV=mV_me*ME_GEV; T=T_me*ME_GEV
    product=G_me_minus2*mV_me*mV_me
    gprime=product/gchi_eff
    gamma_rest=(gprime*gprime)*mV/(12*math.pi)
    factor=analytic_time_dilation(mV,T)
    H=1.66*math.sqrt(gstar)*T*T/MPL_GEV
    return {
     "mV_over_me":mV_me,"T_over_me":T_me,
     "G_times_me2":G_me_minus2,
     "gchi_effective_assumed":gchi_eff,
     "gprime_required":gprime,
     "fixed_coupling_product":product,
     "Gamma_nunu_rest_GeV":gamma_rest,
     "MB_time_dilation_K1_over_K2":factor,
     "Gamma_MB_thermal_per_Zprime_GeV":gamma_rest*factor,
     "H_radiation_reference_GeV":H,
     "Gamma_MB_over_H":gamma_rest*factor/H,
     "existing_MB_vector_Pauli_decay_lower_over_H_at_mu_nu0":gamma_rest*factor/(4*H),
     "conditional_MB_inverse_detailed_balance_not_quantum_FD_BE":True,
     "gstar_assumed_not_computed":gstar,
     "no_thermal_Zprime_abundance_or_physical_neff":True
    }


if __name__=="__main__":
    import json
    print(json.dumps({str(m):conditional_neutrino_decay(mV_me=m)
                      for m in (.03,.1,.3)},indent=2))
