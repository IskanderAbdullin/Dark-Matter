import unittest

import numpy as np
from scipy.special import roots_jacobi

from early_pole_endpoint_warm_positive_gauss_ladder import orthogonal_volume_gauss,volume_gauss


class OrthogonalVolumeGauss(unittest.TestCase):
    def test_high_rules_have_positive_interior_exact_moment_support(self):
        a=np.array([0.,.001,3.,40.])
        b=np.array([1.,.01,3.001,40.00001])
        for n in (7,8):
            q,w,e=volume_gauss(a,b,n)
            self.assertLess(e,3e-13)
            self.assertTrue(np.all((q>a[:,None])&(q<b[:,None])&(w>0.)))
            exact=(b-a)*(a*a+a*b+b*b)/(6*np.pi**2)
            np.testing.assert_allclose(w.sum(axis=1),exact,rtol=2e-14)

    def test_zero_lower_shell_matches_independent_gauss_jacobi(self):
        for n in (7,8):
            z,h=roots_jacobi(n,0.,2.)
            u,w,e=volume_gauss(np.array([0.]),np.array([1.]),n)
            np.testing.assert_allclose(u[0],(1+z)/2,rtol=0,atol=5e-14)
            np.testing.assert_allclose(w[0]/w[0].sum(),h/h.sum(),rtol=0,atol=1e-14)


    def test_old_low_order_abscissae_reproduced_by_different_basis(self):
        a=np.array([0.,.5])
        b=np.array([1.,.7])
        for n in (2,3,4):
            u,w=orthogonal_volume_gauss(a,b,n)
            q,g,_=volume_gauss(a,b,n)
            np.testing.assert_allclose(q,a[:,None]+(b-a)[:,None]*u,rtol=2e-11)
            np.testing.assert_allclose(w,g/g.sum(axis=1)[:,None],rtol=2e-11)


if __name__=='__main__':
    unittest.main()
