import unittest
import numpy as np

from early_pole_endpoint_warm_stress_conditional_hull import hull,values


class StressHullTests(unittest.TestCase):
    def test_massless_linear_integrand(self):
        lo,hi,e,w=hull(1.,2.,0.,0.,.3,16)
        np.testing.assert_allclose(lo,np.full(2,1.3/3))
        np.testing.assert_allclose(hi,lo)
        np.testing.assert_array_equal(e,np.zeros(2))

    def test_off_grid_positive_measure(self):
        q0=np.array([1.03,1.19,1.81])
        p=np.array([.2,.3,.5])
        m0=.07; m=.51
        z=np.hypot([1.,2.],m0)
        t=(p@np.hypot(q0,m0)-z[0])/(z[1]-z[0])
        lo,hi,e,w=hull(1.,2.,m,m0,t,73)
        v=p@values(q0*q0,m)
        self.assertTrue(np.all(lo<=v) and np.all(v<=hi))
        self.assertTrue(np.all(w[0]>=lo) and np.all(w[1]<=hi))


if __name__=='__main__':
    unittest.main()
