"""Analytic two-node tests decoupled from archive-specific collision data."""
import unittest

import numpy as np

from early_pole_daughter_minimum_slope_threshold import solve_joint
from early_pole_daughter_slope_limited_witness import relaxed_outer


class LinearSlopeCheck(unittest.TestCase):
    A=np.array([[.5,.5]])
    m=np.array([[.5],[.5]])
    K=np.array([1.,0.])
    x=np.array([0.,1.])

    def test_independent_minimum(self):
        # Equal cell averages imply both two-node slopes exactly vanish at the
        # feasibility optimum. To change the first occupancies by t requires
        # a shared slope at least t, attained by opposite symmetric tilts.
        for target in (None,.1,.25,.4):
            ans=solve_joint(self.A,self.m,self.K,self.x,None,target)
            self.assertAlmostEqual(ans['verified_witness_slope'],target or 0.,places=11)
            self.assertAlmostEqual(ans['numerical_finite_lower_slope'],target or 0.,places=11)

    def test_weak_dual_box(self):
        # With fixed average 1/2 and cap L, the first occupation is in
        # [1/2-L/2, 1/2+L/2] for L<=1. No chord integration error here.
        for L in (.1,.5,1.):
            lo,hi=relaxed_outer(self.A,self.m[0],self.K,self.x,L,np.array([0.]),0.)
            self.assertAlmostEqual(lo,.5-.5*L,places=11)
            self.assertAlmostEqual(hi,.5+.5*L,places=11)


if __name__=='__main__':
    unittest.main()
