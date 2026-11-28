"""Isolate finite-recoil changes of statistical weights at held spectral commutator."""
import json
from pathlib import Path
import numpy as np

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_spacelike_density_statistics_reference import exact
from finite_mass_qed_spacelike_exact_statistical_aperture_reference import cell
from finite_mass_qed_timelike_pole_batched_absorption_reference import load


def run():
    p=Path(__file__).parent/'results'
    z=json.loads((p/'finite-mass-qed-spacelike-off-root-patch-2026-11-28.json').read_text())['data']
    arch={g:load(g) for g in ('fine','third')}
    out=[]
    for r in z:
        grid,mode,k,w=[r[n] for n in ('grid','mode','k','omega')]
        d,ch,_=arch[grid]
        n0=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
        if mode=='T':
            moments=[cell(C.m,k,w,C.e,C.f) for C in ch]
        else:
            moments=[exact(C.m,C.E,C.f,k,w)[:3] for C in ch]
        G,L,D=np.sum(moments,axis=0)
        if min(G,L,D)<=0 or abs(G-L-D)>3e-14*max(G,L):
            raise ArithmeticError('Exact free directional finite-recoil balance failed')
        A,B,S=r['raw_gain_inverse_signed']
        rho=r['rho']
        retained=r['raw_off_root_estimate_divided_by_n_scalar']
        current=rho*((S+A/n0)*L/D-B)
        direct=rho*(A*(1+1/n0)*L-B*G)/D
        err=abs(direct/current-1)
        if not 0<current or err>3e-12:raise ArithmeticError('Off-root sensitive sign/identity gate')
        out.append(dict(grid=grid,k=k,omega=w,mode=mode,parent_energy_interval=r['parent_energy_interval'],
                        held_commutator_relative_shape_shift=current/retained-1,
                        reconstructed_positive_net_identity_defect=err,
                        finite_free_lesser_over_greater=L/G))
    return dict(scope='Shape-only finite recoil with fixed retarded spectral commutator; not a new dressed propagator.',data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
