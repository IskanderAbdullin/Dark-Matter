"""Test the force inverse over a small positive finite event measure."""
import unittest
import numpy as np
from fd_portal_frozen_response_inverse_reference import solve,potential

class FrozenInverse(unittest.TestCase):
    def test_actual_force_recovery_even_with_correlated_columns(self):
        q=np.array([3.3,3.4,3.5]);g=np.array([.1,.2,.3]);u=np.array([.12,-.15])
        _,F,_=potential(q,g,*u)
        v,n,r=solve(q,g,F,(0.,0.))
        np.testing.assert_allclose(u,v,atol=3e-10,rtol=0)
        self.assertLess(r,2e-13)
