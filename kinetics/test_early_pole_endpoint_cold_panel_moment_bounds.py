import unittest

import numpy as np

from early_pole_endpoint_cold_panel_moment_bounds import bounds


class MonotoneProbes(unittest.TestCase):
    def test_exact_endpoint_bounds_for_native_masses(self):
        p=bounds([0.,1.,3.],[2.,2.],1.)
        np.testing.assert_allclose(p[0],[.5,2.])
        np.testing.assert_allclose(p[1],[np.sqrt(.5),np.sqrt(5.)])
        np.testing.assert_allclose(p[2],[.5/np.sqrt(2.),.5/np.sqrt(2.)+1.5/np.sqrt(10.)])

    def test_zero_total_raises(self):
        with self.assertRaises(ArithmeticError):
            bounds([0.,1.],[0.],1.)


if __name__=='__main__':
    unittest.main()
