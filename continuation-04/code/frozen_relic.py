"""Frozen relativistic BE relic, NR Bose shape, and abundance-convention audit.
Assumptions: one real scalar, mu=0, SM entropy release, no cold condensate.
Run: python frozen_relic.py. Not a Boltzmann perturbation solver.
"""
import json
import numpy as np
from scipy.special import zeta
from scipy.constants import Boltzmann,hbar,c,elementary_charge,pi

def moments():
    m1=3*zeta(4)/zeta(3)
    m2=12*zeta(5)/zeta(3)
    m4=360*zeta(7)/zeta(3)
    return dict(mean_pc_over_T=m1,mean_pc_squared_over_T2=m2,
        mean_pc_fourth_over_T4=m4,
        sigma1D_over_cT_over_mc2=np.sqrt(m2/3),
        kurtosis_excess=(9/5)*m4/(m2*m2)-3)

def abundance(m_eV=86.,gstar=61.75,Tgamma=2.7255,gS_today=3.909):
    ngamma=2*zeta(3)/(pi*pi)*(Boltzmann*Tgamma/(hbar*c))**3/1e6
    nb=(gS_today/gstar)*ngamma/2
    n_source=ngamma/gstar
    rhocrit_h2_eV_cm3=1.87847e-26*c*c/elementary_charge/1e6
    return dict(temperature_CMB_K=Tgamma,one_real_scalar=True,m_eV=m_eV,
      gstarS_dec=gstar,gstarS_today=gS_today,
      n_photon_cm3=ngamma,n_boson_standard_cm3=nb,
      n_boson_source_DMabundance_TcTd1_cm3=n_source,
      ratio_standard_to_source=nb/n_source,
      Omega_h2_standard=m_eV*nb/rhocrit_h2_eV_cm3,
      Omega_h2_source_no_condensate=m_eV*n_source/rhocrit_h2_eV_cm3,
      disclaimer="Different definitions of entropy DOF/condensate factors may explain discrepancy; check original conventions.")

def tails():
    out=[]
    for F in (10,100,1000,10000,100000,1000000):
        L=np.log1p(1/F)
        out.append(dict(F=F,relative_p_volume_aboveF_rel_BE=L**3,
              relative_p_volume_aboveF_NR_critical_BE=L**1.5))
    return out

def run():
    r=dict(frozen_rel_momentum_moments=moments(),reference_abundance=abundance(),
     momentum_only_high_f_volumes=tails(),
     nonrel_free_expansion="If initially BE in p^2, then T~a^-2, fugacity z=constant, n*lambda_th^3=constant; no BEC crossing solely by homogeneous free expansion.",
     non_equivalence="A freeze-out distribution e^(p/T) cannot identically equal a late thermal e^(p^2/(2mT)-mu/T) over a nonzero p interval.",
     caveat="This does not rule out gravitating collapse, nonthermal production, interactions, symmetry charge, or coherent BEC.")
    assert abs(r["frozen_rel_momentum_moments"]["mean_pc_over_T"]-2.7011780329)<1e-9
    assert abs(r["frozen_rel_momentum_moments"]["kurtosis_excess"]-2.0728443082)<1e-9
    assert abs(r["reference_abundance"]["ratio_standard_to_source"]-3.909/2)<1e-12
    return r
if __name__=="__main__":print(json.dumps(run(),indent=2))
