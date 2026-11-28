"""Check that a fixed warm 300-node receipt cannot change its trained chart."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_single_upgrade_history import build


CONTROL=Path(__file__).parent/'results/early-pole-warm-single-g9-upgrade-2026-11-28.json'


class Fixed300(unittest.TestCase):
    def test_chart_is_bound_to_exact_control(self):
        d,m=build(CONTROL)
        self.assertEqual(d.ns,300)
        self.assertEqual(int(np.sum(d.repeats==9)),1)
        self.assertTrue(np.all(d.ds.w>0.))
        with self.assertRaises(ArithmeticError):
            chart_for(m)
        with TemporaryDirectory() as r:
            changed=Path(r)/'control.json'
            changed.write_bytes(CONTROL.read_bytes()+b'\n')
            with self.assertRaises(ArithmeticError):
                chart_for(m,budget_control=changed)
        with self.assertRaises(ArithmeticError):
            chart_for(dict(m,scalar_node_budget=301),budget_control=CONTROL)

    def test_positive_count_and_pressure_envelope(self):
        d,m=build(CONTROL)
        n=np.zeros(d.n)
        n[:d.ns]=1e-8*np.exp(-d.ds.p)
        a=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
        v,b,_=gauss(dict(m,x_end=.0075,steps=4096,nodal_occupation=n),a,.0008,.0075,expected_steps=4096,budget_control=CONTROL)
        self.assertAlmostEqual(v[0]/np.dot(d.ds.w,n[:d.ns]),1,places=14)
        self.assertLessEqual(b['pressure_total_min_max'][0],v[1])
        self.assertGreaterEqual(b['pressure_total_min_max'][1],v[1])
