import unittest

import numpy as np

from early_pole_endpoint_cold_spatial_population_exact_rectangle import decompose,neutrality


class Rectangle(unittest.TestCase):
    def test_separable_total_change_and_fixed_kernel_profile(self):
        raw=np.array([[[2.,3.],[1.,9.]],[[4.,6.],[2.,18.]]])
        r=decompose(raw)
        self.assertAlmostEqual(r['explicit_scalar_share_of_log'],.5)
        self.assertAlmostEqual(r['native_profile_TV'],.3)
        self.assertAlmostEqual(r['symmetric_kernel_profile_TV'],0.)
        self.assertAlmostEqual(r['normalized_profile_interaction_TV'],0.)

    def test_shape_path_interaction_but_additive_signed_contrast(self):
        z=np.array([[[4.,1.],[1.,4.]],[[2.,3.],[2.,3.]]])
        r=decompose(z)
        np.testing.assert_allclose(np.asarray(r['symmetric_scalar_profile_contrast'])+r['symmetric_kernel_profile_contrast'],[-.4,.4])
        self.assertIsNone(r['explicit_scalar_share_of_log'])
        self.assertAlmostEqual(r['normalized_profile_interaction_TV'],.6)

    def test_equal_forward_and_inverse_means_give_equal_roots(self):
        vals=dict(forward_number=6.,inverse_number=2.,forward_excess=3.,inverse_excess=1.)
        fields={k:{'response_rectangle':[[v,v],[v,v]]} for k,v in vals.items()}
        r=neutrality(fields,[1.,2.])
        np.testing.assert_equal(r['number_root'],[[.5,.5],[.5,.5]])
        np.testing.assert_equal(r['energy_root_over_number_root_minus_one'],[[0.,0.],[0.,0.]])



if __name__=='__main__':
    unittest.main()
