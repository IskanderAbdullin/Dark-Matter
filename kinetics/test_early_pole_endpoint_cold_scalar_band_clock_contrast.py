"""Arithmetic fixtures for common and differential late clock biases."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_scalar_band_clock_contrast import contrasts


class ScalarBandClock(unittest.TestCase):
    def test_common_clock_factor_cancels_in_chart_comparison(self):
        u=np.array([2.,3.,5.])
        v=np.array([3.,2.,7.])
        w=np.array([.95,1.04,.98])
        q=contrasts(w*u,u,w*v,v)
        assert_allclose(q['change_in_contrast_fine_minus_coarse'],[0,0,0],rtol=0,atol=4e-16)

    def test_two_distinct_clock_factors_appear_in_chart_ratio(self):
        z=contrasts([1.8],[2.0],[3.3],[3.0])
        self.assertAlmostEqual(z['change_in_contrast_fine_minus_coarse'][0],1.5-(3.3/1.8),places=14)

    def test_rejects_nonpositive_coefficient(self):
        with self.assertRaises(ArithmeticError):
            contrasts([0.],[1.],[2.],[3.])


if __name__=='__main__':
    unittest.main()
