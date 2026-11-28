import unittest
from itertools import product

import numpy as np

from early_pole_endpoint_warm_budgeted_gauss_probe import choose


class BudgetProbe(unittest.TestCase):
    def test_integer_optimizer_agrees_with_small_exhaustive_search(self):
        order=(4,5,7)
        C=np.array([[.2,.1,.09],[.5,.47,.01],[.03,.08,.02]])
        for cap in range(12,22):
            reference=min((sum(C[i,order.index(q)] for i,q in enumerate(p)),sum(p),p)
                           for p in product(order,repeat=3) if sum(p)<=cap)
            n,(loss,path)=choose(C,cap,order)
            self.assertAlmostEqual(loss,reference[0],places=14)
            self.assertEqual((n,path),reference[1:])

    def test_impossible_budget_is_not_silently_raised(self):
        with self.assertRaises(ArithmeticError):
            choose([[.2,.1],[.2,.1]],7,(4,5))


if __name__=='__main__':
    unittest.main()
