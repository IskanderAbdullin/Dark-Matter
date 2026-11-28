"""Sharp source hull, not a pressure or trajectory bound."""
import unittest
import numpy as np
from early_pole_endpoint_warm_actual_source_moment_hull import hull

class SourceMomentHull(unittest.TestCase):
    def test_two_witnesses_at_same_count_and_mean(self):
        u=np.array([0.,1.,2.])
        c=np.array([0.,1.,0.])
        lo,hi=hull(u,c,1.)
        self.assertEqual(lo,(0.,0,2,.5))
        self.assertEqual(hi[0],1.)
        for v,i,j,t in (lo,hi):
            self.assertEqual((1-t)*u[i]+t*u[j],1.)
            self.assertEqual((1-t)*c[i]+t*c[j],v)

    def test_outside_mean_is_rejected(self):
        with self.assertRaises(ArithmeticError):
            hull([0.,1.],[0.,1.],1.1)
