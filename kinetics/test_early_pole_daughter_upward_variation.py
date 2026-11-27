import unittest

import numpy as np

from early_pole_daughter_minimum_upward_variation import solve


class TwoNodeVariation(unittest.TestCase):
    def test_exact_minimum(self):
        A=np.array([[.5,.5]])
        m=np.array([[.3],[.7]])
        K=np.array([1.,0.])
        x=np.array([0.,1.])
        base=solve(A,m,K,x,None,None)
        self.assertAlmostEqual(base['total_upward_variation'],0.,places=10)
        # Native first occupation is at most .6. A +.05 rate thus
        # requires local first occupation <=.55 and local last >=.85.
        pos=solve(A,m,K,x,None,.05)
        self.assertAlmostEqual(pos['total_upward_variation'],.3,places=10)
        self.assertAlmostEqual(pos['numerical_dual_lower'],.3,places=10)


if __name__=='__main__':
    unittest.main()
