"""Thermal-number condition for the sign of the stimulated 1↔2 term."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from early_pole_decay_mode_clock_reference import Q
from early_pole_joint_fermion_thermal_reference import Joint


def one(s):
    x=s['x_end']
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=Joint(shape[tuple(s['cells'])],s['mass_ratio'],0.,'none')
    b=np.asarray(s['endpoint_fermion_occupations'])
    p=np.asarray(s['endpoint_parent_occupations'])
    n_b=d.d.wf@b
    target,pars=d.target(x,b,p)
    theta=-1/pars[2]
    mu=(pars[:2]*theta)
    m=d.d.gap(x,s['endpoint_scalar_occupations'])
    E=np.hypot(m,Q)
    # Common T: fP(Ep)>fD(Ep-Ephi) iff muP-muD>Ephi.
    chem_gate=mu[1]-E
    count_gate=d.d.wf@expit((chem_gate-d.daughter_energies(x))/theta)
    Penergy=np.linspace(.03,8.,121)+E
    diff=lambda cb:expit((mu[1]-Penergy)/theta)-expit((cb-(Penergy-E))/theta)
    assert max(abs(diff(chem_gate)))<3e-15
    assert np.max(diff(chem_gate+.02))<0 and np.min(diff(chem_gate-.02))>0
    sym=expit(-Penergy/theta)-expit(-(Penergy-E)/theta)
    assert np.max(sym)<0
    return dict(x=x,mass_ratio=s['mass_ratio'],joint_temperature_over_T=theta,
                parent_mu_over_T=float(mu[1]),daughter_mu_over_T=float(mu[0]),
                instantaneous_scalar_energy_over_T=E,critical_daughter_mu_over_T=float(chem_gate),
                initial_daughter_number_over_T3=n_b,thermal_daughter_number_at_zero_stimulation_over_T3=count_gate,
                daughter_number_factor_to_stimulus_gate=count_gate/n_b,
                sampled_equal_mu_zero_stimulus_negative=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('specs',nargs='+')
    a=p.parse_args()
    d={}
    for path in a.specs:
        f=Path(path)
        b=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d[f.name]=one(json.loads(b))
    print(json.dumps(d,indent=2,allow_nan=False))
