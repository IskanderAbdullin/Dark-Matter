"""Check averaged log replacement factors and their order sensitivity."""
from itertools import product
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_point_factorial import factors


class LogFactorial(unittest.TestCase):
    def test_symmetric_interaction_is_split_among_factors(self):
        F=np.array([[1.,2.],[.1,.3],[-.4,.2]])
        I=np.array([.2,.4])
        J=np.array([-.6,.9])
        y={k:np.exp(np.array(k)@F+k[0]*k[1]*I+np.prod(k)*J) for k in product((0,1),repeat=3)}
        main,span=factors(y)
        ideal=F+J/3
        ideal[:2]+=I/2
        assert_allclose(main,ideal,atol=1e-15,rtol=0.)
        assert_allclose(main.sum(axis=0),F.sum(axis=0)+I+J,atol=1e-15,rtol=0.)
        self.assertTrue(np.all((span[0]<=main)&(main<=span[1])))

    def test_needs_all_crossed_states(self):
        with self.assertRaises(ArithmeticError):
            factors({(0,0,0):np.ones(2)})


if __name__=='__main__':
    unittest.main()
