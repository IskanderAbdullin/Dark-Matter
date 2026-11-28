import unittest
import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_frozen_ballistic_shape import eval_meas,loss,measure


class BallisticShapeTests(unittest.TestCase):
    def test_massless_single_speed(self):
        for z in (1e-7,.03,.12,1.,2.,9.):
            t,w=roots_legendre(64)
            independent=float(w@np.sin(z*t/2)**2)
            self.assertAlmostEqual(float(loss(z)),independent,places=14)
        y,p,*_=measure(np.array([1.,9.]),np.array([.2,.8]),0.)
        self.assertAlmostEqual(eval_meas(y,p,1.),1.-np.sin(1.))

    def test_equal_second_moment_not_equal_shape(self):
        x=np.array([.5]); p=np.array([1.])
        y=np.array([.2,.8]); q=np.array([.65,.35])
        self.assertAlmostEqual(p@x**2,q@y**2)
        self.assertGreater(q@y**4,p@x**4)
        self.assertNotAlmostEqual(eval_meas(x,p,1.),eval_meas(y,q,1.),places=6)
        self.assertGreater(eval_meas(x,p,1.),eval_meas(y,q,1.))


if __name__=='__main__':
    unittest.main()
