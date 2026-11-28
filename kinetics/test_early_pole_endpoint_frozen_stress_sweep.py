"""Analytic relativistic and deep nonrelativistic moment limits."""
from math import fsum
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_frozen_stress_sweep import evaluate


class FrozenStress(unittest.TestCase):
    def test_massless_equals_rho_over_three_and_pseudo_pressure(self):
        q=np.array([.001,.1,1.,4.])
        w=np.array([2.,8.,.2,.3])
        full,tail=evaluate(q,w,q>=.01,0.)
        assert_allclose(full[:2],[fsum(w),fsum(w*q)],rtol=8e-16,atol=0)
        assert_allclose(full[2:],[full[1]/3,full[1]/3],rtol=8e-16,atol=0)
        self.assertGreater(tail[2]/full[2],tail[0]/full[0])

    def test_nonrelativistic_scaled_stresses(self):
        q=np.array([.001,.1,1.,4.])
        w=np.array([2.,8.,.2,.3])
        m=1e8
        f,_=evaluate(q,w,q>=.01,m)
        assert_allclose(f[1]/m,f[0],rtol=2e-15,atol=0)
        assert_allclose(3*m*f[2],fsum(w*q**2),rtol=2e-15,atol=0)
        assert_allclose(3*m**3*f[3],fsum(w*q**4),rtol=4e-15,atol=0)


if __name__=='__main__':
    unittest.main()
