"""First-slow-coupling response on a fixed archived held elastic orbit."""
import argparse
import gzip
import json
from math import fsum

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
from early_pole_screened_elastic_discrete_gradient_reference import step


def compare(bundle):
    with gzip.open(bundle, 'rt') as stream:
        old = json.load(stream)['branches']
    a, b = [old[k] for k in ('steps', 'cell-conservative')]
    native, continuous = [np.asarray(v['snapshots'][0]['hard_moments']) for v in (a, b)]
    cold = np.asarray(a['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(cold, b['snapshots'][0]['cold_bins'])
    window = a['duration']/a['steps']
    if window != b['duration']/b['steps']:
        raise ValueError('Different first archive windows.')
    C = Screened(a['grid'], 1, len(cold), window, (1, 1, 2), .08, 1.,
                 a['aperture_grid'], 'steps')
    M = C.M
    np.testing.assert_array_equal(native, M.init)
    np.testing.assert_array_equal(cold, C.ns)
    for x,y,z in zip(M.edge_sets,a['hard_edge_sets'],b['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
        np.testing.assert_array_equal(x,z)
    edge = np.asarray(M.edge_sets[2])
    cell = np.flatnonzero((edge[:-1] >= 1.4138888888888888-2e-14)
                          & (edge[1:] <= 1.4666666666666666+2e-14))
    if len(cell) != 2 or not np.allclose(edge[[cell[0], cell[-1]+1]],
                                         (1.4138888888888888, 1.4666666666666666),
                                         rtol=0., atol=2e-14):
        raise ValueError('Require the original exact two daughter cells.')
    local = native.copy()
    ii = M.r[2].start + cell
    local[ii] = continuous[ii]
    h = window/64
    histories, defects = [], []
    for x in (native, local):
        M.init = x.copy()
        elastic = ElasticChannels(C, a['grid'], 1.)
        y = x.copy()
        rates, ledger = [], []
        for j in range(65):
            cold_rhs = Cold.rhs(C, np.r_[y.ravel(), cold, 0., 0.], False)[C.j0:-2]
            rates.append(cold_rhs.tolist())
            if j < 64:
                y, receipt = step(elastic, y, h, 'backward', energy_budget=5e-15)
                ledger.append(receipt['rho_drift'])
        if max(map(abs, ledger)) >= 5e-15:
            raise ValueError('Elastic gate failed on fixed orbit.')
        defects.append(dict(max_abs=max(map(abs, ledger)), signed_sum=fsum(ledger)))
        histories.append(rates)
    delta = np.subtract(histories[1], histories[0])
    report = {}
    for n in (16, 32):
        k = 64//n
        end = delta[::k]
        mid = delta[k//2::k]
        left_end = (end[:-1]+end[1:])/2
        fine = delta[::k//2]
        fine_blocks = ((fine[:-1]+fine[1:])/2).reshape(n, 2, len(cold)).mean(axis=1)
        mixture = (left_end+mid)/2
        residual = float(np.max(np.abs(fine_blocks-mixture)))
        if residual > 8*np.finfo(float).eps*np.max(np.abs(delta)):
            raise ValueError('Composite trapezoid identity failed.')
        report[str(n)] = dict(original_trapezoid_mean=float(left_end.sum()/n),
                              alternate_midpoint_mean=float(mid.sum()/n),
                              mixture_mean=float(mixture.sum()/n),
                              refined_trapezoid_mean=float(fine_blocks.sum()/n),
                              max_block_per_mode_algebra_residual=residual)
    return dict(scope='First derivative in the strength of the entire slow field at zero, on an unchanged 64-leaf held elastic orbit. Not a full-slow extrapolation.',
                leaf_width=h, fixed_seed=cold.tolist(), active_daughter_cells=ii.tolist(),
                elastic_defects=defects, cases=report, cold_source_by_mode=histories)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    print(json.dumps(compare(p.parse_args().bundle), indent=2, allow_nan=False))
