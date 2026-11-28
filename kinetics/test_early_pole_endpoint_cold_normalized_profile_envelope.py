"""Check fractional-bin extrema and certified total variation on finite grids."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_normalized_profile_envelope import compare,shares


class ColdNormalizedProfiles(unittest.TestCase):
    def test_correlated_bin_and_total(self):
        # The more luminous second endpoint puts all stimulated weight in bin 2.
        c=np.array([.2,.2])
        rows=np.array([[.3,0.],[0.,.6]])
        base=c+rows.mean(axis=0)
        p=shares(c,[rows],base)
        assert_allclose(p,[[.2,.5/.7],[.2/.7,.8]],rtol=3e-15,atol=0)

    def test_certificate_counts_only_certain_signs(self):
        a=np.array([[.65,.7],[.3,.35]])
        b=np.array([[.25,.30],[.7,.75]])
        u=dict(count_only=a,count_and_mean_excess=a)
        v=dict(count_only=b,count_and_mean_excess=b)
        q=compare(u,v,np.array([1.,2.,3.]))
        self.assertAlmostEqual(q['count_only']['robust_total_variation_lower_bound'],.35,places=15)
        assert_allclose(q['count_only']['guaranteed_loss_by_common_bin'],[.35,0],atol=3e-16)
        assert_allclose(q['count_only']['guaranteed_gain_by_common_bin'],[0,.35],atol=3e-16)


if __name__=='__main__':
    unittest.main()
