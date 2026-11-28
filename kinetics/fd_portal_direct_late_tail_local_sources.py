"""Resolve residual direct number and heat flow after the equilibrium-clock labels."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_expansion_tracking_continuation import follow
from fd_portal_fast_clock_tail_conserved_coordinates import integrate
from fd_portal_tabulated_entropy_split_expansion import Stage, hubble
from sm_dark_entropy_exchange_reference import StandardModel


def sources(s, p, T0=.90):
    x, T, e, l = p[:4]
    A = np.exp(3*x)
    rho_v = s.eos.state(T)[0]
    rho_d = s.bath([e, l])['rho']
    H = hubble(rho_v, rho_d)
    N0 = s.bath([0., np.log(T0)])['n']
    E0 = s.energy(0., np.array([T0, 0., np.log(T0)]))
    S0 = s.entropy(0., np.array([T0, 0., np.log(T0)]))
    z, sigma = s.collision(T, np.exp(l), e)
    return np.array([-A*z[0]/(N0*H), -A*z[1]/(E0*H), A*sigma/(S0*H)])


def run(eos_path, previous):
    b = previous.read_bytes()
    old = json.loads(b)
    h = hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if old['table_input_sha256'] != h or [r['operator'] for r in old['rows']] != ['timelike','combined']:
        raise ArithmeticError('Need unchanged EOS and both verified direct checkpoints')
    x = np.array([.96, .99, 1.02, 1.05, 1.08, 1.11, 1.14, 1.17, 1.20])
    eos = StandardModel(eos_path)
    pp, qq = PositivePortal((13,36,26)), PositivePortal((17,44,32))
    rows = []
    for record, C in zip(old['rows'], (Stage,CombinedExpansion), strict=True):
        a = np.asarray(record['points'][-1])
        if a[0] != x[0] or a[1] >= record['fast_clock_local_T']:
            raise ArithmeticError('Restart does not lie beyond the previously reported local temperature label')
        s, t = C(eos,pp,112), C(eos,qq,168)
        u, ul = follow(s,a,x,'BDF',1e-11,maximum_step=.0015)
        v, vl = follow(t,a,x,'Radau',2e-12,maximum_step=.0015)
        w = integrate(t,a,x,maximum_step=.0015)
        ext = np.array([[r['log_a'],*r['state'],r['relative_entropy_added'],
                         r['normalized_work'],r['comoving_number_ratio']-1.] for r in w])
        fields = [0,1,2,3,5,6,7]
        d1 = np.max(abs(u-v),axis=0)
        d2 = np.max(abs(v[:,fields]-ext),axis=0)
        if max(d1) > 3e-7 or max(d2) > 4e-8:
            raise ArithmeticError('Independent histories do not resolve the late trajectory')
        z = np.array([sources(t,r) for r in v])
        z2 = np.array([sources(s,r) for r in u])
        if np.min(z[:,2]) < -2e-13 or np.min(np.diff(v[:,5])) < -2e-13 or np.any(np.diff(v[:,1]) >= 0.):
            raise ArithmeticError('Failed entropy or temperature ordering check')
        rows.append(dict(operator=record['operator'],points=v.tolist(),normalized_source_samples=z.tolist(),
                         post_label_delta=(v[-1,[5,6,7]] - v[0,[5,6,7]]).tolist(),
                         source_mesh_defect=np.max(abs(z-z2),axis=0).tolist(),intensive_mesh_defect=d1.tolist(),
                         independent_extensive_defect=d2.tolist(),intensive_ledgers=[ul.tolist(),vl.tolist()],
                         extensive_ledgers=np.max([r['entropy_energy_absolute_defects'] for r in w],axis=0).tolist()))
    return dict(previous_segment_sha256=hashlib.sha256(b).hexdigest(),table_input_sha256=h,fields=old['fields'],
                source_fields=['count_per_log_a','dark_heat_per_log_a','total_entropy_per_log_a'],
                delta_fields=['entropy','pressure_work','count'],rows=rows,
                scope='Direct restricted FD vacuum charged-lepton tail. Earlier local linearized clock temperatures are comparison labels, not decoupling thresholds.')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('previous',type=Path)
    a = p.parse_args()
    print(json.dumps(run(a.eos,a.previous),indent=2,allow_nan=False))
