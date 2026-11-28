import unittest
import numpy as np

from early_pole_endpoint_cold_clock_profile_rectangle import row


class ClockRectangleTests(unittest.TestCase):
    def test_fixed_common_bin_triangle(self):
        q=np.array([1.,2.,3.])
        weights=[[.8,.2],[.7,.3],[.82,.18],[.69,.31]]
        books=[(q,np.array(x).reshape(2,1,1)) for x in weights]
        a=row(q,books,0,0)
        np.testing.assert_allclose(a['S2_and_S4_single_chart_clock_TV'],[.02,.01])
        np.testing.assert_allclose(a['all_cross_TV'],[[.1,.11],[.12,.13]])
        self.assertAlmostEqual(a['maximum_clock_mixing_uncertainty'],.03)


if __name__=='__main__':
    unittest.main()
