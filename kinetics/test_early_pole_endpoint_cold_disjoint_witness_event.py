import unittest

from early_pole_endpoint_cold_disjoint_witness_event import witness


class ExactDisjointWitness(unittest.TestCase):
    def test_nonadjacent_intervals_give_sharp_guaranteed_gap(self):
        q=[0.,1.,2.,3.]
        r=witness(q,q,[.4,.2,.4],[.2,.6,.2],[[0.,1.,0.],[1.,0.,1.]])
        self.assertEqual(r['fine_event_segments'],[[0.,1.],[2.,3.]])
        self.assertEqual(r['third_event_segments'],[[1.,2.]])
        self.assertAlmostEqual(r['fine_event_minimum_fine_minus_third'],.4)
        self.assertAlmostEqual(r['third_event_minimum_third_minus_fine'],.4)

    def test_uncovered_positive_length_intersection_fails(self):
        with self.assertRaises(ArithmeticError):
            witness([0.,1.],[0.,1.],[1.],[1.],[[0.],[0.]])


if __name__=='__main__':
    unittest.main()
