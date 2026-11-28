"""Test disjoint and overlapping normalized profiles against analytic minima."""
import unittest

import numpy as np

from early_pole_endpoint_cold_profile_exact_finite_distance import solve


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


if __name__=='__main__':
    unittest.main()
