"""Number-neutral and photon-heat-neutral temperatures need not be equal."""
import gzip,json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from finite_mass_qed_crossed_annihilation_source_reference import Q,inverse_tail_bound,value
from finite_mass_qed_scalar_absorption_reference import gain_tail_bound


def evaluate(data,mesh):
    measure=[]
    raw=value(data,1.,68.,mesh[0],-1,mesh[1],4,capture=measure)
    k,A,B,_=np.concatenate(measure,axis=0).T
    mass=data['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(Q,mass)
    fs=raw['scalar_f']
    p=float(max(data['endpoint_parent_occupations']))
    b=float(max(data['endpoint_fermion_occupations']))
    if 1-p-b-p*b/fs<=0.:
        raise ArithmeticError('Monotonicity fixture does not hold')
    H=(1+1/fs)*A
    W=H-B
    def integrals(T):
        nb=0. if T==0 else np.exp(-k/T)/(-np.expm1(-k/T))
        z=H+nb*W
        N=z.sum()
        Qgamma=k@z
        return N,Qgamma
    tN=brentq(lambda T:integrals(T)[0],0.,1.,xtol=1e-13)
    tE=brentq(lambda T:integrals(T)[1],0.,1.,xtol=1e-13)
    dN,dE=integrals(tN)[1],integrals(tE)[0]
    q_tail=(68+3)*inverse_tail_bound(data,1.,68.)
    p_tail=(56+3+es)*gain_tail_bound(data,56.)['absolute_gain_tail_upper_bound']
    if max(q_tail,p_tail)>5e-8:raise ArithmeticError('Weighted Planck tail too large.')
    return dict(T_number_neutral=float(tN),T_photon_heat_neutral=float(tE),
                photon_heat_divided_by_held_scalar_occupation_at_T_number_neutral=float(dN),
                scalar_net_divided_by_held_scalar_occupation_at_T_heat_neutral=float(dE),
                charge_crossed_inverse_heat_tail_bound_unscaled=float(q_tail),
                parent_scalar_absorption_gain_heat_tail_bound_unscaled=float(p_tail))


def run(grids):
    tracks=[]
    for g in grids:
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        a,b=(evaluate(d,mesh) for mesh in ((3,4),(6,8)))
        sl=('T_number_neutral','T_photon_heat_neutral')
        drift=max(abs(a[k]-b[k]) for k in sl)
        if drift>5e-6 or abs(b[sl[0]]-b[sl[1]])<.001:
            raise ArithmeticError(('Heat/cardinality root distinction unresolved',g,drift,b))
        tracks.append(dict(grid=g,meshes=[a,b],root_absolute_mesh_drift=drift))
    return dict(scope='Only the IR-finite charged-pair crossing, frozen charges and scalar mode with daughter-shaped anti proxy, solely photon T varied.',data=tracks)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
