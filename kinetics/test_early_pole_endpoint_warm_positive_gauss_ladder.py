"""Independent positive quadrature normalization and limiting-shell checks."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_scalar_positive_cell_quadrature import volume_rule
from early_pole_endpoint_warm_positive_gauss_ladder import volume_gauss


class GaussLadder(unittest.TestCase):
    def test_same_low_order_against_closed_form(self):
        a=np.array([0.,1e-9,.1,4.])
        b=np.array([2e-9,3e-9,.1001,4.8])
        q,w,e=volume_gauss(a,b,2)
        k,v,_=volume_rule(a,b)
        assert_allclose(q,k,rtol=5e-14,atol=0)
        assert_allclose(w,v,rtol=5e-14,atol=0)

    def test_positive_exact_moments_on_scaled_shells(self):
        a=np.array([0.,1e-14,.01,4.])
        b=np.array([2e-10,3e-14,.013,4.8])
        for n in (3,4):
            q,w,_=volume_gauss(a,b,n)
            self.assertTrue(np.all(q>a[:,None]) and np.all(q<b[:,None]) and np.all(w>0))
            for order in range(2*n):
                p=3+order
                expect=(b-a)*sum(a**k*b**(p-1-k) for k in range(p))/(2*np.pi**2*p)
                assert_allclose((w*q**order).sum(axis=1),expect,rtol=5e-14,atol=0)

    def test_higher_volume_moments(self):
        a=np.array([0.,1e-12,.01,4.])
        b=np.array([2e-10,3e-12,.013,4.8])
        q,w,e=volume_gauss(a,b,6)
        self.assertTrue(np.all(w>0) and np.all(q>a[:,None]) and np.all(q<b[:,None]))
        self.assertLess(e,2e-12)
        for k in range(12):
            p=3+k
            expect=(b-a)*sum(a**i*b**(p-1-i) for i in range(p))/(2*np.pi**2*p)
            assert_allclose((w*q**k).sum(axis=1),expect,rtol=2e-12,atol=0)


if __name__=='__main__':
    unittest.main()
