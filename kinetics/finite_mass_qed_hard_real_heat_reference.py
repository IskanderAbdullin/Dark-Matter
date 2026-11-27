"""Separate photon-bath heat from scalar number in two IR-finite real channels."""
import gzip
import json
from pathlib import Path

import numpy as np
from finite_mass_qed_crossed_annihilation_source_reference import value
from finite_mass_qed_scalar_absorption_reference import Q,audit,born


def pair(data,n,angle):
    terms=[]
    z=value(data,1.,68.,n,-1,angle,4,capture=terms)
    k,A,B,C=np.concatenate(terms,axis=0).T
    nb=np.exp(-k)/(-np.expm1(-k))
    fs=z['scalar_f']
    e=np.hypot(Q,data['checkpoints'][-1]['scalar_pole_over_T'])
    prod=A*(1+nb); inverse=B*nb
    observed=prod*(1+fs)-inverse*fs
    scalar=float(e*np.sum(observed))
    photons=float(k@observed)
    charged=-float((k+e)@observed)
    unsigned=float(k@(prod*(1+fs)+inverse*fs))
    norm=abs(scalar)+abs(photons)+abs(charged)+unsigned
    test1=abs((scalar+photons+charged)/norm)
    test2=abs(np.sum(observed)/z['net']-1.)
    if max(test1,test2)>3e-12:
        raise ArithmeticError(('Positive annihilation heat measure failed',test1,test2))
    return dict(net=z['net'],scalar_energy_gain=scalar,photon_bath_gain=photons,
                charged_energy_gain=charged,photon_heat_over_unsigned_turnover=photons/unsigned,
                signed_heat_ledger_scaled_defect=(scalar+photons+charged)/norm,
                photon_weighted_scalar_gain=float(k@prod),photon_weighted_scalar_loss=float(k@inverse))


def run(grids):
    result=[]
    for g in grids:
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        p=[pair(d,n,ang) for n,ang in ((3,4),(6,8))]
        c=[audit(d,1.,56.,*mesh,energy_moments=True) for mesh in ((5,5,18),(8,8,30))]
        thermal=[audit(d,T,38.,6,6,30,equilibrium=True,energy_moments=True) for T in (.4,1.)]
        if max(abs(w['photon_heat_over_unsigned_turnover']) for w in thermal)>3e-13:
            raise ArithmeticError('Thermal heat did not close.')
        keys=('photon_bath_gain','scalar_energy_gain','charged_energy_gain')
        drift=max(abs(b[k]/a[k]-1) for a,b in (p,c) for k in keys)
        B=born(d)
        es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
        result.append(dict(grid=g,pair_annihilation=p[-1],parent_scalar_absorption=c[-1],
                           total_photon_heat_over_Born_scalar_energy_at_g_point005=.005**2*(p[-1]['photon_bath_gain']+c[-1]['photon_bath_gain'])/(B*es),
                           scalar_energy_or_number_correction_at_g_point005=.005**2*(p[-1]['net']+c[-1]['net'])/B,
                           one_charge_born_net=B,max_relative_heat_mesh_drift=drift,thermal_absorption=thermal))
    return dict(scope='Two IR-finite bare real crossings only, one frozen scalar momentum, fixed temperature and stored spectra; no background evolution or full thermal heat.',data=result)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
