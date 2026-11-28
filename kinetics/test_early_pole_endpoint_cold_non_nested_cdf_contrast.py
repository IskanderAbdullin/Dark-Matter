import unittest

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import compare,enclosure


class PanelCDF(unittest.TestCase):
    def test_no_interpolation_and_exact_cuts(self):
        z=enclosure(np.array([0.,1.,2.]),[3.,1.],[0.,.4,1.,1.3,2.])
        np.testing.assert_allclose(z,[[0.,0.],[0.,.75],[.75,.75],[.75,1.],[1.,1.]])

    def test_distinct_boundaries_can_certify_a_gap(self):
        r=compare([0.,.4,1.],[0.,.6,1.],[.8,.2],[.1,.9])
        self.assertAlmostEqual(r['normalized_TV_lower_bound_via_cdf'],.7,places=12)
        self.assertEqual(r['strongest_q'],.4)

    def test_no_claim_for_overlapping_enclosures(self):
        r=compare([0.,.4,1.],[0.,.6,1.],[.4,.6],[.6,.4])
        self.assertEqual(r['normalized_TV_lower_bound_via_cdf'],0.)

    def test_sharp_overlap_combines_nonadjacent_discrepancies(self):
        edges=[0.,1.,2.,3.]
        r=compare(edges,edges,[.4,.2,.4],[.2,.6,.2])
        self.assertAlmostEqual(r['normalized_TV_lower_bound_via_cdf'],.2,places=12)
        a,b=r['sharp_compatible_TV']['minimum_TV_bounds']
        self.assertLessEqual(a,.4)
        self.assertLessEqual(.4,b)
        self.assertLess(b-a,1e-12)

    def test_definite_crossing_panel_may_allow_identical_measures(self):
        a=[0.,.25,1.]
        b=[0.,.75,1.]
        r=compare(a,b,[.3,.7],[.6,.4])
        l,h=r['sharp_compatible_TV']['minimum_TV_bounds']
        self.assertEqual(l,0.)
        self.assertLess(h,1e-12)



if __name__=='__main__':
    unittest.main()
