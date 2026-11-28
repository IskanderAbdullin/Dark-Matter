import unittest

import numpy as np

from finite_mass_qed_native_scalar_bridge_join import field,ledger


class PastBridgeLedger(unittest.TestCase):
    def test_separate_conditional_roots_preserve_total_energy(self):
        n=np.array([[3.,1.,2.],[1.5,.2,1.3],[1.2,.1,1.1]])
        w=.8
        h=field(n,w,.1)
        self.assertNotAlmostEqual(h[0,1]/h[0,2],h[1,1]/h[1,2])
        self.assertLess(abs(h[2,0]-(w*n[0,0]-n[1,0])),1e-14)
        r=ledger([n],w,.1)
        self.assertTrue(all(abs(sum(t['scalar_photon_and_charged_energy']))<1e-14 for t in r['tests']))

    def test_inconsistent_energy_refuses_join(self):
        n=np.array([[3.,1.,2.],[1.5,.4,1.1],[1.2,.1,1.1]])
        with self.assertRaises(ArithmeticError):
            field(n,.8,.1)


if __name__=='__main__':
    unittest.main()
