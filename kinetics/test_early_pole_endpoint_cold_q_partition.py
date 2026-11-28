"""Synthetic exact-grid checks for normalization and conditional redistribution."""
import unittest

import numpy as np

from early_pole_endpoint_cold_q_partition import neutral,regroup,shape


class ColdQPartition(unittest.TestCase):
    def test_coarsened_profile_difference_does_not_need_interpolation(self):
        e=np.array([1.,2.,3.])
        s=shape(e,np.array([1.,3.]),e,np.array([3.,1.]),e)
        self.assertEqual(s['TV_lower_bound_from_exact_common_bins'],.5)
        self.assertEqual(s['Wasserstein1_q_panel_supported_interval'][0],0.)
        self.assertAlmostEqual(sum(x['signed_normalized_difference'] for x in s['important_bins']),0.)

    def test_count_neutral_total_preserves_opposite_bin_fluxes(self):
        e=np.array([1.,2.,3.])
        z=np.zeros((2,3,3))
        z[:,0,1]=[1.,3.]
        z[:,0,2]=[2.,2.]
        r=neutral(e,z,e)
        self.assertEqual(r['number_neutral_occupation'],1.)
        self.assertEqual(r['positive_flux_sum_over_total_gross_inverse'],.25)
        self.assertEqual(r['negative_flux_sum_over_total_gross_inverse'],-.25)
        self.assertEqual(r['sum_roundoff'],0.)

    def test_no_nonshared_boundary(self):
        with self.assertRaises(ArithmeticError):
            regroup(np.array([1.,2.,3.]),np.ones(2),np.array([1.,2.5,3.]))


if __name__=='__main__':
    unittest.main()
