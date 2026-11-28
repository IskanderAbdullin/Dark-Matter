"""Cross-check minimax witnesses and bounds by tiny exhaustive enumerations."""
from itertools import product
import unittest

import numpy as np
from early_pole_endpoint_warm_budget_minimax import fit, values


class Minimax(unittest.TestCase):
    def test_double_epoch_gap(self):
        loss=np.array([[[1.,0.],[0.,.1]],[[0.,.1],[1.,0.]]])
        o=(4,5)
        cand,v,(lo,hi),w=fit(loss,o,9)
        self.assertLessEqual(sum(cand),9)
        self.assertAlmostEqual(max(v),1.)
        self.assertLessEqual(lo,1.)
        self.assertAlmostEqual(hi,1.)
        self.assertAlmostEqual(sum(w),1.)
        self.assertLess(lo,hi)

    def test_matches_enumeration(self):
        rng=np.random.default_rng(1423)
        orders=(4,5,6)
        feasible=[p for p in product(orders,repeat=3) if sum(p)<=15]
        for _ in range(8):
            q=rng.uniform(0.,1.,(4,3,3))
            p,v,(lo,hi),w=fit(q,orders,15)
            best=min(max(values(q,orders,path)) for path in feasible)
            self.assertLessEqual(lo,best+1e-14)
            self.assertAlmostEqual(hi,best)
            self.assertAlmostEqual(max(v),best)
            self.assertIn(p,feasible)
            self.assertAlmostEqual(sum(w),1.)


if __name__=='__main__':
    unittest.main()
