"""The q-panel support enclosure must not interpolate any subpanel data."""
import unittest
import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_band_panel_contrast import bound


class PanelEnclosure(unittest.TestCase):
    def test_straddling_cut_keeps_all_unknown_weight(self):
        edges=np.array([1e-6,4e-6,5e-6,1e-5,1.5e-5])
        q,p=bound(edges,np.array([10.,20.,60.,10.]))
        assert_allclose(q,[[1e-6,4e-6],[5e-6,1e-5],[5e-6,1e-5],[1e-5,1.5e-5]],rtol=0,atol=1e-20)
        assert_allclose(p,[[.9,1.],[.7,.7],[.1,.7],[.1,.1],[0,.1]],rtol=0,atol=1e-15)
    def test_negative_or_overlap_rejected(self):
        with self.assertRaises(ArithmeticError):
            bound(np.array([1e-6,4e-6,3e-6]),np.array([1.,2.]))
        with self.assertRaises(ArithmeticError):
            bound(np.array([1e-6,3e-6,4e-6]),np.array([1.,-1.]))


if __name__=='__main__':
    unittest.main()
