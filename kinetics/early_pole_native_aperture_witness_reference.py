"""Exact frozen bare 1↔2 aperture for explicit native moment-matched profiles."""
import argparse
import gzip
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals

GRID = Path(__file__).resolve().parent / 'results'
PROFILES = GRID / 'early-pole-native-aperture-witness-profiles-2026-11-27.json'


def occup(e, nodes, values):
    j = np.searchsorted(nodes, e, side='right') - 1
    return np.where((j >= 0) & (j < len(values)), values[np.clip(j, 0, len(values)-1)], 0.)


def refine(edges, values, cells, mass):
    points = np.unique(np.r_[edges, np.concatenate([c['qEdges'][1:-1] for c in cells])])
    mid = (points[:-1] + points[1:]) / 2
    base = occup(mid, edges, values)
    new = base.copy()
    errors = []
    z, w = roots_legendre(16)
    for c in cells:
        p = np.asarray(c['qEdges'])
        j = np.argmin(np.abs(edges[:-1] - p[0]))
        assert np.allclose(edges[j:j+2], p[[0,-1]], rtol=0., atol=2e-15)
        assert abs(values[j]-c['baseline']) < 2e-15
        f = np.asarray(c['f']).copy()
        assert 0. <= f.min() <= f.max() <= 1.
        mask = (p[0] < mid) & (mid < p[-1])
        q = (p[:-1,None]+p[1:,None])/2 + np.diff(p)[:,None]*z/2
        wq = np.diff(p)[:,None]*w/2 * q*q
        E = np.hypot(q, mass)
        # Undo rounding of the printed profile using three strictly interior
        # occupations. The centered basis spans 1,E,1/E, but avoids their
        # near-degeneracy on a very narrow native cell.
        E0 = np.hypot((p[0]+p[-1])/2, mass)
        dE = (q-p[0])*(q+p[0])/(E+np.hypot(p[0],mass))
        dE0 = ((p[0]+p[-1])/2-p[0])*((p[0]+p[-1])/2+p[0])/(E0+np.hypot(p[0],mass))
        dE -= dE0
        width = np.diff(np.hypot(p[[0,-1]],mass))[0]
        mat = np.array([wq.sum(axis=1), (wq*dE/width).sum(axis=1),
                        (wq*(dE/width)**2*E0/E).sum(axis=1)])
        eligible = np.flatnonzero((f > 1e-3)&(f < 1-1e-3))
        i = max(combinations(eligible,3), key=lambda i:abs(np.linalg.det(mat[:,i])))
        correction = np.linalg.solve(mat[:,i], -mat@(f-values[j]))
        assert np.max(np.abs(correction)) < 1e-9
        f[list(i)] += correction
        assert 0. <= f.min() <= f.max() <= 1.
        new[mask] = occup(mid[mask], p, f)
        dc = f[:,None] - values[j]
        defect = [float(np.sum(wq*dc*measure)) for measure in (1., E, 1/E)]
        assert max(map(abs,defect)) < 4e-18
        errors.append(dict(defects=defect, largest_profile_roundoff_correction=float(np.max(abs(correction)))))
    return points, base, new, errors


def mode(M, m, b, q, p_edges, original, witness, d_edges, d_f, alpha):
    E, lo, hi, _, _ = intervals(M, m, b, q)
    ep = np.hypot(p_edges, M)
    ed = E + np.hypot(d_edges, b)
    mesh = np.unique(np.r_[lo, hi, ep[(lo < ep)&(ep < hi)], ed[(lo < ed)&(ed < hi)]])
    y = (mesh[:-1] + mesh[1:]) / 2
    w = np.diff(mesh)
    P0 = occup(y, ep, original)
    P1 = occup(y, ep, witness)
    D = occup(y, ed, d_f)
    P = P0 + alpha*(P1-P0)
    # Independent conversion into parent radial integration, segmented at all
    # parent+daughter edges. Energy derivative holds for the bare pole only.
    pu = np.sqrt(np.maximum(0., (mesh-M)*(mesh+M)))
    z, wz = roots_legendre(8)
    r = (pu[:-1,None]+pu[1:,None])/2 + np.diff(pu)[:,None]*z/2
    wr = np.diff(pu)/2 * ((r/np.hypot(r,M))@wz)
    delta = max(abs((w-wr)@a) for a in (P,D,P*(1-D)))
    return E, np.array([w@P, w@D, w@(P*(1-D))]), float(delta)


def evaluate(profile):
    label = profile['grid']
    path = GRID / f'early-pole-spectrum-{label}-m06-2026-11-26.json.gz'
    with gzip.open(path, 'rt') as fh:
        a = json.load(fh)
    M = a['x_end']
    m = a['checkpoints'][-1]['scalar_pole_over_T']
    b = a['mass_ratio']*M
    p_edges, original, witness, errors = refine(
        np.asarray(a['endpoint_parent_momentum_edges']),
        np.asarray(a['endpoint_parent_occupations']), profile['parent'], M)
    d_edges = np.asarray(a['endpoint_fermion_momentum_edges'])
    d_f = np.asarray(a['endpoint_fermion_occupations'])
    qs = a['scalar_momentum_nodes']
    v = np.diff(np.asarray(a['scalar_momentum_edges'])**3)/(6*np.pi**2)
    n = np.asarray(a['endpoint_scalar_signed_shell_number_over_T3'])
    alphas = [0., .25, .5, .75, .85, 1.]
    rows = []
    worst = 0.
    for alpha in alphas:
        mode_integrals = [mode(M,m,b,q,p_edges,original,witness,d_edges,d_f,alpha) for q in qs]
        worst = max(worst, *(r[2] for r in mode_integrals))
        k = int(np.argmax(n))
        Es = np.array([r[0] for r in mode_integrals])
        P, D, seed = np.array([r[1] for r in mode_integrals]).T
        net = (n*(P-D)+v*seed)/(Es*qs)
        rows.append(dict(witness_fraction=alpha, net=float(net.sum()), cold_net=float(net[np.array(qs)<1e-5].sum()),
                         seed=float((v*seed/(Es*qs)).sum()), dominant_node=dict(q=qs[k], parent_integral=P[k],
                              daughter_integral=D[k], stimulated=float(P[k]-D[k])),
                         accumulated_parent_change=float((n*P/(Es*qs)).sum())))
    assert worst < 1e-13
    assert rows[0]['net'] > 0. > rows[-2]['net']
    crossover = rows[0]['net']/(rows[0]['net']-rows[-1]['net'])
    for r in rows:
        np.testing.assert_allclose(r['net'], rows[0]['net']*(1-r['witness_fraction']/crossover),
                                   rtol=3e-12, atol=3e-9)
    vol=np.diff(p_edges**3)/3
    near=np.maximum(0., (np.minimum(p_edges[1:],1.48)**3-np.maximum(p_edges[:-1],1.39)**3)/3)
    l1 = [float(z@abs(witness-original)/(z@original)) for z in (vol,near)]
    return dict(grid=label,parent_mass=M,scalar_mass=m,daughter_mass=b,moment_defects_N_E_J=errors,
                full_witness_parent_relative_L1_full_and_near_aperture=l1,
                max_independent_radial_vs_energy_integral_absolute_defect=worst,
                linear_zero_witness_fraction=crossover,rows=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles',type=Path,default=PROFILES)
    args = p.parse_args()
    print(json.dumps(dict(scope='Frozen bare 1↔2; finite archived supports and fixed scalar/daughter. '
                                'Dimensionless kernel omits a common positive prefactor.',
                          cases=[evaluate(c) for c in json.loads(args.profiles.read_text())]),
                     indent=2,allow_nan=False))
