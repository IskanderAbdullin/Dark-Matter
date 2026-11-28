import unittest
import numpy as np

from early_pole_endpoint_cold_forward_total_anchor import solve


def card(v):
    r=np.ones((3,1))
    return np.zeros_like(r),r,[np.asarray(v)],[np.asarray(v)]


class ForwardAnchorTests(unittest.TestCase):
    def test_matching_forward_forces_other_rate_to_move(self):
        a=card([[0.,0.,0.],[2.,1.,0.]])
        b=card([[1.,0.,0.]])
        z=solve(a,b,(0,1,2),1)
        self.assertAlmostEqual(z['minimum_finite_loss'],.5)
        self.assertAlmostEqual(z['support_lower_bound'],.5)
        self.assertAlmostEqual(z['projected_upper_witness'],.5)
        self.assertEqual(z['matched_forward_totals_in_S2_scale'],[1.,1.])

    def test_detects_unreachable_match(self):
        with self.assertRaises(ArithmeticError):
            solve(card([[0.,0.,0.],[.8,1.,0.]]),card([[1.,0.,0.]]),(0,),0)


if __name__=='__main__':
    unittest.main()
