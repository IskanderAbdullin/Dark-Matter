"""Check exact interval-wise L1 allocation, including sign crossings."""
import unittest
from early_pole_endpoint_warm_actual_source_rank_budget import loss

class RankBudget(unittest.TestCase):
    def test_common_window_has_independently_attainable_bounds(self):
        self.assertEqual(loss([[-2.,1.],[2.,4.],[-3.,-2.]]),[4.,9.])
    def test_reversed_endpoint_is_invalid(self):
        with self.assertRaises(ArithmeticError): loss([[2.,1.]])
