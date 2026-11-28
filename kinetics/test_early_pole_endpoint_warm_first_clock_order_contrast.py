import unittest

import numpy as np

from early_pole_endpoint_warm_first_clock_order_contrast import run,values


def card(k,x):
    a=np.asarray([[2.,3.,2.],[2.,3.,2.]])*x
    replay=dict(x_start=.0075,x_end=.01,steps=[512,1024],warm_gauss_order=k,
                fixed_pressure_mass=.2,maximum_trace_absolute_difference=0.,
                coarse_count_and_pressures=a[0].tolist(),fine_count_and_pressures=a[1].tolist(),
                coarse_over_fine_minus_one=[0.]*3)
    compare=dict(x=.01,warm_gauss_order=k,fixed_mass=.2,warm_gauss_moments=a[1].tolist(),S4_reference=[2.,3.,2.])
    return replay,compare


class FirstClock(unittest.TestCase):
    def test_order_only_contrast_has_no_clock_effect(self):
        a,b=card(5,1.),card(6,1.2)
        r=run([a[0],b[0]],[a[1],b[1]])
        np.testing.assert_allclose(r['G6_over_G5_minus_one_at_each_clock'],[[.2]*3]*2)
        np.testing.assert_equal(r['fine_minus_coarse_of_G6_over_G5'],[0.]*3)

    def test_fine_endpoint_must_bind_independent_comparison(self):
        a,b=card(6,1.)
        b['warm_gauss_moments'][0]*=.99
        with self.assertRaises(ArithmeticError):
            values(a,b,6)


if __name__=='__main__':
    unittest.main()
