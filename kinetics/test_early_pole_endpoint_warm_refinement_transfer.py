import unittest
import numpy as np

from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_warm_targeted_source_refinement import PARENTS,Targeted
from early_pole_endpoint_warm_refinement_transfer_test import stats


class TransferRefinementTests(unittest.TestCase):
    def test_exact_boundary_and_repeated_lift(self):
        b=load('fine')[0]
        g=Targeted(b,list(PARENTS))
        a=WarmGauss(b,4)
        self.assertEqual((g.ns,a.ns),(219,211))
        self.assertEqual(g.initial.shape,(g.n,))
        self.assertTrue(np.array_equal(g.initial[g.ns:],a.initial[a.ns:]))
        with self.assertRaises(ArithmeticError):
            Targeted(b,[PARENTS[0],PARENTS[0]])
        with self.assertRaises(ArithmeticError):
            Targeted(b,[PARENTS[0]+1e-8])
        with self.assertRaises(ArithmeticError):
            Targeted(b,[float(a.ds.ed[-1])])

    def test_score_preserves_negative_error(self):
        q={'base':np.array([1.,3.]),'split':np.array([1.8,2.3]),'reference':np.array([2.,2.]),
           'scale':4.,'edges':np.array([1.,2.,3.]),'uniform5':np.array([2.1,2.1])}
        z=stats(q,[0],'selected')
        self.assertAlmostEqual(z['signed'],.2)
        self.assertAlmostEqual(z['sum_abs'],.3)
        self.assertEqual(z['refined_parent_edges'],[1.])


if __name__=='__main__':
    unittest.main()
