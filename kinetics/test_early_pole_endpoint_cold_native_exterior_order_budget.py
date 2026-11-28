"""Test the sharp equal-cap condition independently."""
import unittest
import numpy as np
from early_pole_endpoint_cold_native_exterior_order_budget import fractions

class ExteriorBudget(unittest.TestCase):
    def test_exact_opposite_endpoint_witness(self):
        a=np.array([2.,1.8,1.2])
        b=np.array([1.,1.7,1.1])
        H=np.array([5.,6.,2.])
        e=fractions(a,b,H)
        for q, sign in ((.999,-1),(1.001,1)):
            np.testing.assert_array_equal(np.sign((1-q*e)*(b-a)+q*e*H),np.full(3,sign))

    def test_order_required(self):
        with self.assertRaises(ArithmeticError):
            fractions([1.],[2.],[6.])
