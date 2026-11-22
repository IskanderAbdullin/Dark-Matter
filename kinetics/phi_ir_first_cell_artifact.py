"""Thermal real-boson IR first-cell collocation artifact, not BEC.

The project's finite-volume φ mode measure is exact at cell level,
BUT sampling thermal b(E) at E0=mφ and multiplying by the entire
positive first-cell phase-space capacity spuriously overestimates
Nφ as mφ/h->0. The physical Bose p²dp integral is finite at mφ=0.

This module only diagnoses initial equilibrium measure; the
first-cell averaged Bose occupancy is NOT inserted into the
on-shell Yukawa quantum-collision interpolation (which would need an
independent detailed-balance/H proof). No 86-eV relic prediction.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import zeta
from kinetics.yukawa_bose_pauli_bounded_time import density_cells


def thermal_number_integral(m,T,pmax,g=1):
    """Physical (g/(2π²))∫p²dp BE(sqrt(p²+m²)/T), μ=0.

    Analytic p=0,m=0 integrand is a removable zero (limit 0).
    """
    if not (m>=0 and T>0 and pmax>=0 and g>0 and
            all(math.isfinite(x) for x in (m,T,pmax,g))):
        raise ValueError('finite nonnegative mass/momentum and positive T,g')
    def integrand(p):
        e=math.hypot(p,m)
        x=e/T
        if x==0 or x>700:return 0.
        return p*p/math.expm1(x)
    val,err=quad(integrand,0,pmax,epsabs=2e-13,epsrel=2e-11,limit=250)
    return g/(2*math.pi**2)*val


def thermal_first_cell_and_total(m,*,upper=12.1,bins=64,T=1.,g=1):
    if not(m>0 and upper>m and isinstance(bins,int) and bins>=4 and T>0 and g>0):
        raise ValueError('point-collocated BE requires strictly positive m and finite domain')
    grid=np.linspace(m,upper,bins+1)
    w=density_cells(grid,m,g)
    b=1/np.expm1(grid/T)
    n_node=float(np.dot(w,b))
    h=grid[1]-grid[0]
    emax_first=m+h/2
    pmax_first=math.sqrt(emax_first**2-m*m)
    first_true=thermal_number_integral(m,T,pmax_first,g)
    first_node=float(w[0]*b[0])
    total_true=thermal_number_integral(m,T,math.sqrt(upper*upper-m*m),g)
    corrected=n_node-first_node+first_true
    return {'mass_GeV_if_input_GeV':m,'upper':upper,'temperature':T,
        'bins':bins,'spacing':h,'first_cell_capacity':float(w[0]),
        'first_cell_point_sampled_number':first_node,
        'first_cell_exact_thermal_number':first_true,
        'first_cell_false_enhancement_ratio':first_node/first_true,
        'full_point_sampled_total':n_node,
        'full_exact_thermal_total':total_true,
        'full_false_enhancement_ratio':n_node/total_true,
        'total_if_ONLY_first_cell_number_corrected':corrected,
        'only_first_cell_corrected_ratio':corrected/total_true,
        'first_cell_true_average_occupation':first_true/float(w[0]),
        'asymptotic_first_cell_ratio_h_over_3m':h/(3*m),
        'not_a_real_phi_BEC_or_cosmology':True}


def audit(*,full=False):
    cases=[]
    masses=(.1,.01,.001,.0001,86e-9)
    bins=(16,64,256) if full else (16,64)
    for m in masses:
        for n in bins:
            cases.append(thermal_first_cell_and_total(m,upper=12.1,bins=n))
    massless_n_inf=float(zeta(3)/(math.pi**2))
    return {'schema_version':1,
      'derived_asymptotic':'m<<h<<T: first-bin point sample / physical integral = h/(3*m) [1+O(h/T,m/h)]',
      'massless_BE_finite_number_T1_g1_infinite_cutoff':massless_n_inf,
      'massless_BE_finite_number_T1_g1_cut50':thermal_number_integral(0,1,50),
      'cases':cases,
      'warning':'86eV=86e-9 GeV here, but T=1 GeV and upper=12.1 GeV are illustrative, NOT self-consistent FRW. Large low-mode occupation multiplied by finite first cell volume is NOT BEC.',
      'class_ready':False}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--full',action='store_true')
    x=p.parse_args()
    print(json.dumps(audit(full=x.full),indent=2))
