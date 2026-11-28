"""Check symmetric common-bath and signed non-common FD derivatives directly."""
import unittest
import numpy as np
from fd_portal_differential_reciprocity_reference import PositivePortal,integrate

class DifferentialReciprocity(unittest.TestCase):
    def test_common_state_matches_original_equilibrium_weights(self):
        p=PositivePortal((9,24,16))
        T=.07
        f,D,a,W=integrate(p,T,0.,0.)
        np.testing.assert_array_equal(f,np.zeros(3))
        self.assertEqual(a,0.)
        q=p.moments(T,T,0.)
        np.testing.assert_allclose(D,[[q[3],q[4]/T],[q[4]/T,q[5]/T**2]],rtol=2e-13,atol=0.)
        self.assertGreater(min(W),0.)

    def test_one_signed_affinity_controls_differential_asymmetry(self):
        p=PositivePortal((9,24,16))
        for k in (-1,1):
            f,D,a,_=integrate(p,.07,k*.4,k*.08)
            self.assertGreater(k*a,0.)
            self.assertGreater(f[2],0.)
            self.assertAlmostEqual((D[0,1]-D[1,0])/a,1.,places=12)
