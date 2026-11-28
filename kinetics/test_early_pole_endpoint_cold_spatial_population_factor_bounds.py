import unittest

import numpy as np

from early_pole_endpoint_cold_spatial_population_factor_bounds import crossed,normalized_bounds,tv_with_interval_profile


class CrossPopulation(unittest.TestCase):
    def test_target_cell_edge_is_not_counted_twice(self):
        q=np.array([0.,.5,1.])
        b,m=crossed(q,q,np.array([5.,10.]),q,[4.,9.],q,[1.,3.],False)
        np.testing.assert_allclose(b,[[2.,2.],[4.,4.]])
        self.assertEqual(m,0)

    def test_straddling_population_does_not_interpolate_kernel(self):
        source=np.array([0.,1.])
        target=np.array([0.,.5,1.])
        b,m=crossed(source,source,np.array([6.]),source,[2.],target,[1.,5.],True)
        np.testing.assert_allclose(b,[[3.,15.]])
        self.assertEqual(m,1)
        for x in np.linspace(0.,3.,10):
            self.assertLessEqual(b[0,0],1*x+5*(3-x))
            self.assertLessEqual(1*x+5*(3-x),b[0,1])

    def test_intervals_contain_all_positive_normalized_corners(self):
        B=normalized_bounds(np.array([[2.,4.],[3.,5.]]))
        expected=(2./7.,4./7.)
        np.testing.assert_allclose(B[0],expected)
        for x in np.linspace(2.,4.,4):
            for y in np.linspace(3.,5.,4):
                for s in (.05,.35,.65,.95):
                    p=np.array([s,1-s])
                    lo,hi=tv_with_interval_profile(p,B)
                    a=np.array([x,y])/(x+y)
                    self.assertLessEqual(lo,abs(p-a).sum()/2)
                    self.assertLessEqual(abs(p-a).sum()/2,hi)

    def test_refuse_zero_inverse_denominator(self):
        q=np.array([0.,.5,1.])
        with self.assertRaises(ArithmeticError):
            crossed(q,q,np.ones(2),q,[4.,0.],q,[1.,3.],True)


if __name__=='__main__':
    unittest.main()
