"""Decisive examples for quantile-bin-based Wasserstein support bounds."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_peer_charge_shape_transport import quantile_support_bounds


class PeerSupport(unittest.TestCase):
    def test_touching_support_needs_no_strict_common_atom(self):
        e=np.array([0.,1.,2.])
        a=np.array([0.,.8,1.])
        b=np.array([0.,.2,1.])
        lo,hi,atom=quantile_support_bounds((e,a),(e,b))
        self.assertEqual(lo,0.)
        self.assertIsNone(atom)
        self.assertGreater(hi,1.)

    def test_positive_gap_for_separated_quantile_support(self):
        e=np.array([0.,.7,1.3,2.])
        a=np.array([0.,.8,.9,1.])
        b=np.array([0.,.1,.2,1.])
        lo,hi,atom=quantile_support_bounds((e,a),(e,b))
        assert_allclose(lo,.36,rtol=3e-15)
        self.assertGreater(hi,lo)

    def test_exact_witness_on_overlapping_support(self):
        e=np.array([0.,.7,2.])
        f=np.array([0.,1.,2.])
        p=np.array([0.,.6,1.])
        lo,hi,atom=quantile_support_bounds((e,p),(f,p))
        self.assertEqual(lo,0.)
        self.assertEqual(atom['atoms'],2)
        self.assertEqual(atom['max_absolute_cell_mass_defect'],0.)


if __name__=='__main__':
    unittest.main()
