"""Verify symmetry, zero-count convention and effect isolation."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_warm_pressure_cells import partition


class WarmPressureCells(unittest.TestCase):
    def test_effect_isolation_and_sign_reversal(self):
        a=np.array([[1.,4.],[2.,6.],[1.,3.],[0.,0.],[0.,0.]])
        b=np.array([[2.,8.],[2.,4.],[0.,0.],[1.,3.],[0.,0.]])
        result=np.array([[4.,4.,0.],[-2.,0.,-2.],[-3.,-3.,0.],[3.,3.,0.],[0.,0.,0.]])
        assert_allclose(partition(a,b),result,atol=0,rtol=0)
        assert_allclose(partition(b,a),-result,atol=0,rtol=0)

    def test_general_count_shape_interaction(self):
        a=np.array([[1.,4.]])
        b=np.array([[3.,6.]])
        assert_allclose(partition(a,b),[[2.,6.,-4.]],atol=0,rtol=0)

    def test_empty_cell_cannot_carry_pressure(self):
        with self.assertRaises(ValueError):
            partition(np.array([[0.,1.]]),np.array([[1.,1.]]))


if __name__=='__main__':
    unittest.main()
