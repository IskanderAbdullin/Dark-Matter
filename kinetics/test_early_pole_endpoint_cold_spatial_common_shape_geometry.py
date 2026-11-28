"""Common-shape CDF bounds are not independent panel placements."""
import unittest
import numpy as np
from early_pole_endpoint_cold_spatial_common_shape_geometry import band

class Geometry(unittest.TestCase):
    def test_outer_bounds_are_sharp_on_shared_bin_shape(self):
        e=np.array([0.,1.,2.,4.])
        w=np.array([.1,.2,-.3])
        r=band(e,w)
        np.testing.assert_allclose(r['flow_cdf_L1_q_lower_upper'],[.1,1.])
        self.assertAlmostEqual(r['fixed_common_shape_TV'],.3)

    def test_crossing_allows_zero_contribution(self):
        r=band([0.,1.,2.,3.],[.1,-.2,.1])
        np.testing.assert_allclose(r['flow_cdf_L1_q_lower_upper'],[0.,.3])
        with self.assertRaises(ArithmeticError):
            band([0.,1.],[.1])
