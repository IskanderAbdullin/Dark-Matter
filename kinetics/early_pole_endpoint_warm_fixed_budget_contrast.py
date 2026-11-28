"""Compare actual mixed255 and uniform positive Gauss warm prefixes on one fixed S4 mass."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_budgeted_warm_history import audit
from early_pole_endpoint_gauss_warm_chain_receipt import run as chain, EPOCHS, CLOCK
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_warm_pressure_cells import partition
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart, SPECS, load


def run(control, native_receipt, native_states, mixed_reports, mixed_states, g7_reports, g7_states, *, attribution=False, uniform_order=7):
    if uniform_order not in (7,8):
        raise ArithmeticError("Require a validated uniform order, G7 or G8")
    k = len(native_states)
    if not 1 <= k <= len(EPOCHS) or any(len(seq) != k for seq in
            (mixed_reports, mixed_states, g7_reports, g7_states)):
        raise ArithmeticError('Require the same actual chronological prefix on all three charts')
    a = audit(control, mixed_reports, mixed_states)
    b = chain(g7_reports, g7_states, prefix=True, order=uniform_order)
    out = []
    d = ScalarSubshellChart(load('fine')[0], SPECS['fine'], 2, 4, 1e-5)
    for i, (x, clock) in enumerate(zip(EPOCHS[:k], CLOCK[:k], strict=True)):
        c, h = read(native_states[i])
        q = native_receipt['validated_chain'][i]
        if (h, q['x'], c['x_end'], c['steps']) != (q['uncompressed_state_sha256'], x, x, 2*clock):
            raise ArithmeticError('The common S4 mass does not bind its accepted state')
        occ = np.asarray(c['endpoint_occupation'])
        if occ.shape != (d.n,) or not np.all(np.isfinite(occ)) or min(occ) < 0 or max(occ[d.ns:]) > 1:
            raise ArithmeticError('Invalid native S4 occupations')
        m = d.state(x, occ)[0]
        vals, bounds, cells = [], [], []
        for f, ar, hk, fixed in ((mixed_states[i], a['states'][i], 'state_sha256', control),
                (g7_states[i], b['checked_epoch_pairs'][i], 'uncompressed_state_sha256', None)):
            s, h = read(f)
            if h != ar[hk]:
                raise ArithmeticError('State bytes changed after the audit')
            v, e, cgroup = gauss(s, d, m, x, expected_steps=2*clock, budget_control=fixed)
            vals.append(v)
            cells.append(cgroup)
            bounds.append(e['pressure_total_min_max'])
        l, u = bounds
        p = [l[0]/u[1]-1., l[1]/u[0]-1.]
        contrast = vals[0]/vals[1]-1.
        if not p[0]-2e-14 <= contrast[1] <= p[1]+2e-14:
            raise ArithmeticError('Fixed budget pressure escaped its conditional bounds')
        details={}
        if attribution:
            tab=partition(cells[1],cells[0])/vals[1][1]
            tot=[fsum(row) for row in tab.T]
            if abs(tot[0]-contrast[1])>2e-14 or abs(tot[0]-fsum(tot[1:]))>2e-14:
                raise ArithmeticError('The all-cell pressure partition failed')
            ix=np.argsort(-abs(tab[:,0]))[:5]
            details=dict(cell_partition=dict(fields=['total','changed_count','changed_pressure_per_particle'],
                all_cells=tot,leading=[dict(q=d.ds.ed[j:j+2].tolist(),contributions=tab[j].tolist()) for j in ix]))
        out.append(dict(x=x, steps=[clock, 2*clock], common_native_sha256=q['uncompressed_state_sha256'],
                        fixed_native_mass=m, moments=[v.tolist() for v in vals],
                        **{f"mixed255_over_G{uniform_order}_minus_one":contrast.tolist()}, conditional_pressure_contrast=p,**details))
    return dict(charts=['mixed255',f'G{uniform_order}'], fields=['count','pressure','warm_pressure'],
                budget_reference_sha256=a['budget_reference_sha256'], contrasts=out,
                mixed_trace=a['states'], uniform_trace=b['checked_epoch_pairs'],
                scope='Actual paired prefixes at common S4 measurement mass. Cell count/energy bounds do not include solver error; held-source ranking and physical pressure are not inferred.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ('control','native-chain'):
        p.add_argument('--'+flag, type=Path, required=True)
    for flag in ('native-states','mixed-reports','mixed-states','g7-reports','g7-states'):
        p.add_argument('--'+flag, type=Path, nargs='+', required=True)
    p.add_argument('--uniform-order',type=int,choices=(7,8),default=7)
    p.add_argument('--attribution',action='store_true')
    a = p.parse_args()
    print(json.dumps(run(a.control,json.loads(a.native_chain.read_text()),a.native_states,a.mixed_reports,
                         a.mixed_states,a.g7_reports,a.g7_states,attribution=a.attribution,uniform_order=a.uniform_order),indent=2,allow_nan=False))
