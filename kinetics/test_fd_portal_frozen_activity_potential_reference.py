"""Test the fixed-state test-force potential without the physical state derivative."""
import unittest
import numpy as np
from fd_portal_frozen_activity_potential_reference import potential

class FrozenActivity(unittest.TestCase):
    def test_zero_drive_and_even_convex_potential(self):
        q=np.array([.5,2.,8.])
        g=np.array([.7,.32,.6])
        y,f,H=potential(q,g,0.,0.)
        self.assertEqual(y,0.)
        np.testing.assert_array_equal(f,np.zeros(2))
        self.assertGreater(min(np.linalg.eigvalsh(H)),0.)
        a,b,h=potential(q,g,.3,-.2)
        c,d,k=potential(q,g,-.3,.2)
        self.assertEqual(a,c)
        np.testing.assert_array_equal(b,-d)
        np.testing.assert_array_equal(h,k)
