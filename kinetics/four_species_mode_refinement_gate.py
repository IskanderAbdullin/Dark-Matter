"""Independent Fermi-logit energy mesh error order in coupled χ1χ2φ+lepton.

Only an instantaneous on-shell finite-vector 2↔2 Pauli weak-to-strong
source at fixed incoming 5D quadrature; never infer physical Γelastic.
"""
from __future__ import annotations
import argparse
import json
import numpy as np
from kinetics.four_species_cutoff_mesh_separation import compare_mesh


def linear_interpolation_bound(h, amplitude, wavenumber):
    """Exact interpolation remainder sup|L-Lh| <= h² sup|L''| / 8."""
    if h<=0 or amplitude<0 or wavenumber<=0:raise ValueError
    return amplitude*wavenumber*wavenumber*h*h/8


def evaluate_rate_refinement(bins=(256,512,1024,2048,4096)):
    x=compare_mesh(cut=12.,mesh_cut=20.,bins=bins,
                   nlocal=3,nmu=6,nc=12,nphi=12)
    rows=x['projected_5D_mode_occupations']
    for r in rows:
        r['chi1_logit_uniform_error_bound']=linear_interpolation_bound(
            r['energy_step'],.3,4.)
    entropy_errors=[abs(r['entropy_error_relative']) for r in rows]
    ratios=[entropy_errors[i]/entropy_errors[i+1]
            for i in range(len(entropy_errors)-1)]
    return {'schema_version':1,'stage':'four species fermionic logit mode interpolation order',
            'quad_scope':'same incoming 12x3 Gauss, nmu6, CM mapped12, azimuth12, physical incoming cutoff12, mesh cutoff20; 2 separate analytic-f weak species',
            'reference':x['reference_5D_analytic_occupations'],
            'mesh_rows':rows,
            'H_error_ratios_on_bin_doubling':ratios,
            'theoretical_logit_remainder_bound':'h² sup|logit(f) second derivative|/8',
            'not_a_continuum_5D_or_physical_dm_gamma':True,'class_ready':False}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--full',action='store_true')
    a=p.parse_args()
    grid=(256,512,1024,2048,4096) if a.full else (256,512,1024)
    print(json.dumps(evaluate_rate_refinement(grid),indent=2))
