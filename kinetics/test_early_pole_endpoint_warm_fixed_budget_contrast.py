"""Exercise the chart binding without fabricating a paired trajectory."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from early_pole_endpoint_budgeted_warm_history import build
from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart, SPECS, load

CONTROL = Path(__file__).parent/'results/early-pole-warm-budgeted-gauss-probe-2026-11-28.json'


class TestFixedBudgetContrast(unittest.TestCase):
    def test_requires_frozen_selection_and_exact_byte_digest(self):
        d, card = build(CONTROL)
        with self.assertRaisesRegex(ArithmeticError, 'unbound'):
            chart_for(card)
        label, q, repeats = chart_for(card, budget_control=CONTROL)
        self.assertEqual((label, q.ns, q.nf, q.np), ('mixed255',255,448,488))
        np.testing.assert_array_equal(repeats, d.repeats)
        with TemporaryDirectory() as temp:
            changed = Path(temp)/'control.json'
            changed.write_bytes(CONTROL.read_bytes()+b'\n')
            with self.assertRaisesRegex(ArithmeticError, 'frozen selection'):
                chart_for(card, budget_control=changed)
        broken = dict(card, warm_orders_by_cell=card['warm_orders_by_cell'][::-1])
        with self.assertRaisesRegex(ArithmeticError, 'frozen selection'):
            chart_for(broken, budget_control=CONTROL)

    def test_positive_fixed_geometry_keeps_count_and_pressure_envelope(self):
        d, card = build(CONTROL)
        b = ScalarSubshellChart(load('fine')[0], SPECS['fine'], 2, 4, 1e-5)
        y = np.zeros(d.n)
        y[:d.ns] = 1e-8*np.exp(-d.ds.p)
        z, e, cells = gauss(dict(card,x_end=.0075,steps=4096,nodal_occupation=y),
                             b,.0008,.0075,expected_steps=4096,budget_control=CONTROL)
        self.assertAlmostEqual(z[0]/np.dot(d.ds.w,y[:d.ns]),1.,places=14)
        self.assertAlmostEqual(z[1]/sum(cells[:,1]),1.,places=14)
        lo, hi = e['pressure_total_min_max']
        self.assertLessEqual(lo,z[1])
        self.assertGreaterEqual(hi,z[1])
