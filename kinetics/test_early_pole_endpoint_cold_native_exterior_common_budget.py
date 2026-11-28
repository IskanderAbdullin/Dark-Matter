"""Check asymmetry and shared cancellation with explicit exterior measures."""
import unittest
import numpy as np
from early_pole_endpoint_cold_native_exterior_common_budget import unmatched

class SharedBudget(unittest.TestCase):
    def test_extreme_witness_exact_threshold(self):
        a,b,H=np.array([2.]),np.array([1.]),np.array([5.])
        for lam in (0.,.4,.9):
            w=unmatched(a,b,H,lam)[0]
            for A,B in ((w[0],0.),(0.,w[1]),(w[2],w[2])):
                for shared_probe in (0.,.1,4.):
                    p=(1-lam-A)*a+lam*shared_probe
                    q=(1-lam-B)*b+lam*shared_probe+B*H
                    self.assertLess(abs(p-q)[0],2e-15)
                    self.assertGreater((q+1e-6*B*(H-b)-p+1e-6*A*a)[0],0.)
