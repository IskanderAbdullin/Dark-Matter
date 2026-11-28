import unittest

import numpy as np

from early_pole_endpoint_cold_two_segment_clock_rectangle import decompose


class TwoSegmentClock(unittest.TestCase):
    def test_multiplicative_clock_effects_have_no_log_interaction(self):
        r=decompose([[1.,.9],[.8,.72]])
        self.assertAlmostEqual(r['earlier_only_log_shift'],np.log(.8))
        self.assertAlmostEqual(r['final_only_log_shift'],np.log(.9))
        self.assertAlmostEqual(r['interaction_log_shift'],0.)

    def test_nonseparability_is_retained_per_response(self):
        r=decompose([[[2.,3.],[2.,3.]],[[2.,3.],[4.,1.]]])
        np.testing.assert_allclose(r['interaction_log_shift'],[np.log(2.),np.log(1./3.)])
        with self.assertRaises(ArithmeticError):
            decompose([[0.,.9],[.8,.72]])


if __name__=='__main__':
    unittest.main()
