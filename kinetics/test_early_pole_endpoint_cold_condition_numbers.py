"""Analytical logarithmic perturbation controls independent of the collision callback."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_condition_numbers import estimate


class LocalConditionNumbers(unittest.TestCase):
    def test_power_law_slopes(self):
        f=lambda m,w: np.array([m**2*w**(-4),m**(-3)*w**7])
        z=estimate((1.7,3.1),f)
        assert_allclose(z[0]['richardson_slopes'],[2.,-3.],rtol=4e-11,atol=0)
        assert_allclose(z[1]['richardson_slopes'],[-4.,7.],rtol=4e-11,atol=0)
        assert_allclose(z[0]['local_log_parameter_change_for_one_percent_coefficient'],
                        np.log1p(.01)/np.array([2,3]),rtol=4e-11,atol=0)


if __name__=='__main__':
    unittest.main()
