"""Ensure child projection preserves original panel mass without hiding topology."""
import unittest
import numpy as np
from early_pole_endpoint_cold_native_placement_refinement import metrics, project


class NativePlacement(unittest.TestCase):
    def test_project_and_bound(self):
        e=np.array([0.,1.,2.])
        f=np.array([0.,.5,1.,2.])
        base=np.array([.2,.8])
        a=np.array([.05,.17,.81])
        b,tv,drift=project(e,base,f,a)
        np.testing.assert_allclose([sum(b[:2]),b[2]],base,atol=1e-15)
        self.assertAlmostEqual(b[0]/b[1],a[0]/a[1])
        self.assertGreater(tv,0.)
        self.assertAlmostEqual(drift,.1)
        m0=metrics(f,f,a,b,1.)
        self.assertLessEqual(m0['minimum_q_bracket'][1],2*tv+1e-12)

    def test_reject_crossed_child(self):
        with self.assertRaises(ArithmeticError):
            project(np.array([0.,1.,2.]),np.ones(2),np.array([0.,.5,2.]),np.ones(2))


if __name__=='__main__':
    unittest.main()
