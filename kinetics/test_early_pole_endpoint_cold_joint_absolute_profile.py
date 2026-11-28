import unittest
import numpy as np

from early_pole_endpoint_cold_joint_absolute_profile import solve


def card(vertices, offset=None):
    origin = np.ones((3, 1))
    return (np.zeros_like(origin) if offset is None else np.array(offset)[:, None],
            origin, [np.array(vertices)], [np.array(vertices)])


class JointAbsoluteResponseTest(unittest.TestCase):
    def test_individually_reachable_is_not_jointly_reachable(self):
        a = card([[0., 1., .5], [1., 0., .5]])
        b = card([[0., 0., .5]])
        self.assertAlmostEqual(solve(a, b, (0,), 0)['minimum_l1'], 0.)
        self.assertAlmostEqual(solve(a, b, (1,), 0)['minimum_l1'], 0.)
        self.assertAlmostEqual(solve(a, b, (0, 1, 2), 0)['minimum_l1'], 1.)

    def test_spontaneous_offset_cannot_be_renormalized_away(self):
        a = card([[2., 1., 2.]], [1., 0., 1.])
        b = card([[2., 1., 2.]])
        r = solve(a, b, (0, 1, 2), 1)
        self.assertEqual(r['witness_l1_by_field'], [1., 0., 1.])
        self.assertAlmostEqual(r['support_hyperplane_value'], 2.)
        self.assertAlmostEqual(r['projected_witness_value'], 2.)
        self.assertLess(r['max_verification_defect'], 1e-12)

    def test_independent_bound_mode(self):
        a=card([[0.,1.,.5],[1.,0.,.5]])
        b=card([[0.,0.,.5]])
        r=solve(a,b,(0,1,2),0,method='highs-ipm',independent_bracket=True)
        self.assertAlmostEqual(r['support_hyperplane_value'],1.)
        self.assertAlmostEqual(r['projected_witness_value'],1.)
        self.assertLess(r['verified_bracket_width'],1e-9)


if __name__ == '__main__':
    unittest.main()
