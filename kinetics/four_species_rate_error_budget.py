"""Independent axes of the finite χ1χ2 Pauli 5D rates in the four-pop network.

Use the CURRENT FD modes and moving real-φ Yukawa source, but all rates
are synthetic dimensionless toy parameters. This is a numerical error
BUDGET, not physical Γelastic/H or a validated 86-eV relic density.
"""
from __future__ import annotations
import argparse
import json
import math
from kinetics.four_species_yukawa_shared_lepton import configuration,initial,rhs


def sample(*,cut=12.,nlocal=3,nmu=8,nc=16,nphi=16,mlepton=0.,
           bins=32,T1=1.1,T2=.9,Tlep=1.,Tphi=1.,
           mesh_cut=None,shape='FD'):
    """Orthogonalize the *incoming collision cutoff* from common energy mesh.

    For mesh_cut=None the historical two limits move together; the comparison
    is a documented negative control. mesh_cut>=cut fixes logit grid and
    Yukawa mode capacities while changing only incoming 5D phase space.
    """
    segments=int(math.ceil(cut))
    c=configuration(bins=bins,elastic_cut=cut if mesh_cut is None else mesh_cut,
       elastic_integration_cut=cut if mesh_cut is not None else None,
       mlepton=mlepton,
       elastic_radial_scheme='segmented_p',elastic_segments=segments,
       elastic_nr=nlocal,elastic_nmu=nmu,elastic_nc=nc,elastic_nphi=nphi)
    x=initial(c,T1=T1,T2=T2,Tlep=Tlep,Tphi=Tphi)
    if shape=='oscillatory':
        import numpy as np
        from scipy.special import expit
        g1,g2=c['elastic']['chi_grid'];lg=c['elastic']['lg']
        x=(expit(-g1/T1+.3*np.sin(4*g1)),
           expit(-g2/T2-.2*np.sin(3*g2)),
           expit(-lg/Tlep+.25*np.cos(4*lg)),x[3])
    elif shape!='FD':
        raise ValueError('unsupported synthetic occupation shape')
    _,o=rhs(c,x)
    return {'cutoff_T':cut,'mesh_cutoff_T':c['elastic']['cut'],
       'shape':shape,'energy_segments':segments,'local_radial_order':nlocal,
       'incoming_relative_angle':nmu,'CM_log_transfer_order':nc,'azimuth':nphi,
       'mlepton':mlepton,'bins':bins,'Tchi1':T1,'Tchi2':T2,
       'energy_grid_hchi1':float(c['elastic']['chi_grid'][0][1]-c['elastic']['chi_grid'][0][0]),
       'Qchi1_elastic':o['elastic_Qchi'][0],
       'Qchi2_elastic':o['elastic_Qchi'][1],
       'Qlepton_elastic':o['elastic_Qlep'],
       'entropy_elastic':o['elastic_entropy'],
       'entropy_yukawa':o['yukawa_entropy'],
       'combined_entropy_H_node_error':o['entropy_identity_abs_error'],
       'combined_energy_residual':o['energy_residual'],
       'max_charge_residual':max(abs(z) for z in o['charge_residuals']),
       'Yukawa_Rphi':o['Yukawa_Rphi'],
       'no_physical_DM_86eV_Gamma_or_CLASS':True}


def audit(axis, *, full=False):
    choices={
      'cutoff':('cut',(10.,12.,16.,20.) if full else (10.,12.,16.)),
      'incoming_mu':('nmu',(4,6,8,10,12,16) if full else (4,6,8)),
      'outgoing_cm':('nc',(8,12,16,24) if full else (8,12,16)),
      'outgoing_phi':('nphi',(8,12,16,24) if full else (8,12,16)),
      'radial_local':('nlocal',(2,3,4,5) if full else (2,3,4)),
      'lepton_mass':('mlepton',(0.,.5)),
    }
    if axis not in choices:raise ValueError('unknown independent rate axis')
    param,values=choices[axis]
    results=[]
    for value in values:
        kw={param:value}
        # Compare one numerical axis at fixed all other inputs.
        results.append(sample(**kw))
    return {'schema_version':1,'axis':axis,'rows':results,
            'scope':'synthetic χ1χ2 full-Pauli 5D with moving real φ Yukawa; rates not dimensionful physical Γel, no SM/FRW/86eV relic',
            'class_ready':False}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--axis',default='cutoff',choices=('cutoff','incoming_mu','outgoing_cm','outgoing_phi','radial_local','lepton_mass'))
    p.add_argument('--full',action='store_true')
    a=p.parse_args()
    print(json.dumps(audit(a.axis,full=a.full),indent=2))
