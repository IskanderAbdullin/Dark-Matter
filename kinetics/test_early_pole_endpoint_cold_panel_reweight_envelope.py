"""Exact synthetic checks on the moment polytope and correlated neutrality."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_panel_reweight_envelope import envelope,mean_vertices


class ColdPanelReweights(unittest.TestCase):
    def test_exact_mean_vertices(self):
        x=np.array([1.,2.,3.])
        a=mean_vertices(x,2.)
        assert_allclose(a,[[.5,0.,.5],[0.,1.,0.]],atol=0,rtol=0)
        assert_allclose(a@x,[2.,2.],atol=0,rtol=0)

    def test_joint_neutrality_not_independent_division(self):
        c=np.array([1.,0.,1.])
        options=np.array([[2.,1.,1.],[6.,3.,3.]])
        full=c+options.mean(axis=0)
        v=envelope(c,[options],full)
        assert_allclose(v['conditional_neutrality_relative_to_original_min_max'],[.75,1.125],rtol=2e-14)
        assert_allclose(v['coefficient_relative_to_original_min_max'],
                        [[.6,.5,2/3],[1.4,1.5,4/3]],rtol=2e-14)


if __name__=='__main__':
    unittest.main()
