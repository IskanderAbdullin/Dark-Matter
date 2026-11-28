"""Test disjoint and overlapping normalized profiles against analytic minima."""
from itertools import product
import unittest

import numpy as np

from early_pole_endpoint_cold_profile_exact_finite_distance import (
    normalized_support_lower, profile_from_projected_blocks, solve,
)


def fixed(p):
    return dict(offset=np.array(p),base=np.array(p),total=1.,count_only=[])


def interval():
    c=np.array([.2,.2])
    r=np.array([[.3,0.],[0.,.6]])
    bas=c+r.mean(axis=0)
    return dict(offset=c,base=bas/bas.sum(),total=bas.sum(),count_only=[r])


class JointProfileLP(unittest.TestCase):
    def test_disjoint_profiles(self):
        r=solve(fixed([.9,.1]),interval(),'count_only')
        self.assertAlmostEqual(r['minimum_finite_common_bin_TV'],.9-(.5/.7),places=10)

    def test_can_match_exactly_in_overlap(self):
        r=solve(fixed([.6,.4]),interval(),'count_only')
        self.assertAlmostEqual(r['minimum_finite_common_bin_TV'],0.,places=10)
        self.assertAlmostEqual(r['direct_witness_TV'],0.,places=10)

    def test_inverse_style_zero_offset_retains_uniformly_positive_total(self):
        c=dict(offset=np.zeros(2),base=np.array([.35,.65]),total=1.,
               count_only=[np.array([[.3,.7],[.4,.6]])])
        r=solve(fixed([.8,.2]),c,'count_only')
        self.assertAlmostEqual(r['minimum_finite_common_bin_TV'],.4,places=10)

    def test_independent_bracket_matches_analytic_separation_and_overlap(self):
        for A,expected in (([.9,.1],.9-.5/.7),([.6,.4],0.),([.1,.9],.2-.1)):
            for a,b in ((fixed(A),interval()),(interval(),fixed(A))):
                r=solve(a,b,'count_only',independent_bracket=True)['independent_profile_bracket']
                self.assertLessEqual(r['support_lower_bound'],expected+1e-15)
                self.assertLessEqual(expected,r['projected_profile_upper_bound'])
                self.assertLess(r['verified_TV_bracket_width'],4e-13)

    def test_lifted_support_matches_enumerated_multicell_vertices(self):
        off=np.array([1.2,0.,.35])
        blocks=[np.array([[.1,.4,.9],[.6,0.,.4]]),
                np.array([[.5,.2,0.],[0.,.2,1.4],[0.,.1,.35]])]
        card=dict(offset=off,total=2.,count_only=blocks)
        for w in (np.array([-.5,.425,.05]),np.array([0.,0.,0.]),np.array([.2,-.3,-.5])):
            raw=[off+a+b for a,b in product(*blocks)]
            exact=min(w @ (p/p.sum()) for p in raw)
            lower=normalized_support_lower(card,'count_only',w)
            self.assertLessEqual(lower,exact+1e-15)
            self.assertLess(exact-lower,4e-13)

    def test_empty_and_negative_optimizer_blocks_still_give_a_witness(self):
        c=interval()
        q=profile_from_projected_blocks(c,'count_only',np.array([0.,-.3,0.]))
        np.testing.assert_allclose(q,c['base'],rtol=0.,atol=2e-16)



if __name__=='__main__':
    unittest.main()
