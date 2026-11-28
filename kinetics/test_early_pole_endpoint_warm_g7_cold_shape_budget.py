"""Enclose adversarial cold releases with a direct extremal construction."""
import unittest

import numpy as np

from early_pole_endpoint_warm_g7_cold_shape_budget import release_loss


class ColdShapeBudget(unittest.TestCase):
    def test_common_to_independent_loss_is_exact(self):
        left=np.array([.4,.9,.1,.2])[:,None]
        right=np.array([.8,.1,.1,.7])[:,None]
        A=np.array([[1.,2.],[.3,.6],[.8,.8],[.2,1.]])
        B=np.array([[2.,6.],[2.,4.],[2.,5.],[3.,3.]])
        shared=np.minimum((right-left)*A,(right-left)*B)
        independent=right*A-left*B
        np.testing.assert_allclose(shared-independent,release_loss(left,right,A,B),rtol=4e-15)
        for eps in (.03,.4,1.):
            rhs=shared-eps*release_loss(left,right,A,B)
            lhs=(1-eps)*shared+eps*independent
            np.testing.assert_allclose(lhs,rhs,rtol=4e-15)
