import unittest
import numpy as np

from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_warm_selective_order_control import SelectiveOrder


class SelectiveOrderTests(unittest.TestCase):
    def test_paired_cells_only(self):
        cfg=load('fine')[0]
        g=WarmGauss(cfg,4)
        for k,n in ((5,213),(6,215)):
            a=SelectiveOrder(cfg,k)
            self.assertEqual((a.ns,a.nf,a.np),(n,448,488))
            self.assertEqual(sum(a.repeats!=g.repeats),2)
            self.assertEqual(a.initial.shape,(a.n,))
            ia=np.r_[0,np.cumsum(a.repeats)]
            ig=np.r_[0,np.cumsum(g.repeats)]
            for j in np.flatnonzero(a.repeats==g.repeats):
                np.testing.assert_array_equal(a.ds.p[ia[j]:ia[j+1]],g.ds.p[ig[j]:ig[j+1]])
                np.testing.assert_array_equal(a.ds.w[ia[j]:ia[j+1]],g.ds.w[ig[j]:ig[j+1]])


if __name__=='__main__':
    unittest.main()
