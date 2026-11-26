"""Orthogonal cutoff and Fermi logit energy-grid numerical audit.

Both χ1 and χ2 exchange energy with ONE leptonic Pauli reservoir and
real φ still evolves under moving χ1<->χ2+φ Yukawa event when time is
advanced. Here we audit its INSTANTANEOUS full 5D elastic weak moments.

The historical scan changed the integration cutoff and the discrete
energy grid simultaneously. For FD, logit(f) is affine in energy, so
changing grid size alone does not bias the reconstructed occupation.
For NONTHERMAL f, grid interpolation can alias strongly. Keep common
mesh fixed while changing incoming phase-space cutoff and compare
against a separate direct 5D evaluation of the analytic spectrum.

This is a numerical diagnostic with arbitrary units/couplings, NOT
cosmological Γelastic, ΩDM or a 86eV CLASS-ready boson spectrum.
"""
from __future__ import annotations
import argparse
import json
import numpy as np
from scipy.special import expit
from kinetics.four_species_rate_error_budget import sample
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def analytic_occupation(M, *, T1=1.1, T2=.9, Tlep=1.):
    if M not in (.1, .02):
        raise ValueError('only two declared illustrative Dirac χ masses')
    def f1(E):
        E=np.asarray(E)
        return expit(-E/T1+.3*np.sin(4*E))
    def f2(E):
        E=np.asarray(E)
        return expit(-E/T2-.2*np.sin(3*E))
    def fl(E):
        E=np.asarray(E)
        return expit(-E/Tlep+.25*np.cos(4*E))
    return (f1 if M==.1 else f2),fl


def independent_weak(cut=12.,nlocal=3,nmu=6,nc=12,nphi=12):
    rows=[]
    for M in (.1,.02):
        f,g=analytic_occupation(M)
        x=evaluate((f,g),M=M,m=0.,g=.01,mv=.1,
            nr=nlocal,nmu=nmu,nc=nc,nphi=nphi,cutoff=cut,
            radial_scheme='segmented_p',segments=int(cut),
            angular_scheme='log_transfer')
        rows.append({'mass':M,'Qchi':x['Qchi'],'entropy':x['entropy'],
            'energy_defect':x['max_energy_defect'],
            'negative_entropy_events':x['negative_entropy_events']})
    return rows


def compare_mesh(*,cut=12.,mesh_cut=20.,bins=(32,64,128,256),
                 nlocal=3,nmu=6,nc=12,nphi=12):
    """Projection error on SAME physical quadrature, NO rate calibration."""
    weak=independent_weak(cut,nlocal,nmu,nc,nphi)
    source=[]
    for b in bins:
        s=sample(cut=cut,mesh_cut=mesh_cut,bins=b,shape='oscillatory',
                 nlocal=nlocal,nmu=nmu,nc=nc,nphi=nphi)
        source.append({'bins':b,'energy_step':s['energy_grid_hchi1'],
           'Qchi1':s['Qchi1_elastic'],'Qchi2':s['Qchi2_elastic'],
           'Qlepton':s['Qlepton_elastic'],
           'entropy':s['entropy_elastic'],
           'Qchi1_error_relative':s['Qchi1_elastic']/weak[0]['Qchi']-1,
           'Qchi2_error_relative':s['Qchi2_elastic']/weak[1]['Qchi']-1,
           'entropy_error_relative':s['entropy_elastic']/sum(x['entropy'] for x in weak)-1,
           'energy_residual':s['combined_energy_residual'],
           'entropy_identity_error':s['combined_entropy_H_node_error']})
    return {'mesh_cut':mesh_cut,'incoming_cut':cut,
            'incoming_radial':f'{int(cut)} segments x {nlocal}',
            'incoming_mu':nmu,'CM_log_angle':nc,'azimuth':nphi,
            'reference_5D_analytic_occupations':weak,
            'projected_5D_mode_occupations':source,
            'scope':'synthetic prescribed nonthermal Fermi mode alias under full on-shell finite-mediator Pauli with ONE shared lepton; no cosmological Γelastic',
            'CLASS_ready':False}


def audit(full=False):
    cut_grid=[8.,12.,16.] if not full else [8.,10.,12.,16.,20.]
    records=[]
    for cut in cut_grid:
        a=sample(cut=cut,bins=64,shape='oscillatory',nlocal=3,nmu=6,nc=12,nphi=12)
        b=sample(cut=cut,mesh_cut=20.,bins=64,shape='oscillatory',nlocal=3,nmu=6,nc=12,nphi=12)
        records.append({'incoming_cut':cut,'moving_mesh_Qchi1':a['Qchi1_elastic'],
            'fixed_mesh_Qchi1':b['Qchi1_elastic'],
            'moving_mesh_H':a['entropy_elastic'],
            'fixed_mesh_H':b['entropy_elastic'],
            'relative_Qchi1_bias_due_to_mesh':a['Qchi1_elastic']/b['Qchi1_elastic']-1,
            'relative_H_bias_due_to_mesh':a['entropy_elastic']/b['entropy_elastic']-1})
    out=compare_mesh(bins=(32,64,128,256) if full else (32,64,128))
    return {'schema_version':1,'mesh_and_cutoff_confounded_scan':records,
            'grid_refinement_at_same_incoming_cutoff':out,
            'warning':'H/N/E perfect source may still underestimate full nonthermal continuum-event entropy by tens of percent at coarse energy grid',
            'physical_86eV_relic_validated':False}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--full',action='store_true')
    a=p.parse_args()
    print(json.dumps(audit(full=a.full),indent=2))
