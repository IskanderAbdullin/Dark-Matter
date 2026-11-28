import unittest
import numpy as np

from early_pole_endpoint_warm_bounded_shape_witness import bands


class BoundedShapeTests(unittest.TestCase):
    def test_same_count_and_energy_in_both_families(self):
        for mean in (.0001,.21,.5,.8,.9999):
            for width in (.05,.10,.20):
                for split in (False,True):
                    z,w=bands(mean,width,split)
                    self.assertTrue(np.all(z>=0) and np.all(z<=1))
                    self.assertTrue(np.all(np.diff(z,axis=1)>0))
                    self.assertAlmostEqual(w.sum(),1.)
                    self.assertAlmostEqual(w@z.mean(axis=1),mean)


if __name__=='__main__':
    unittest.main()
