"""Independent interval-support transport fixtures."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_cold_profile_cdf_support import distances


class CdfPanelBounds(unittest.TestCase):
    def test_touching_panels_do_not_imply_nonzero_distance(self):
        e=np.array([0.,1.,2.])
        a=dict(base=[1.,0.])
        b=dict(base=[0.,1.])
        r=distances(e,e,a,b,'fixed')
        assert_allclose(r['Wasserstein1_q_lower_upper'],[0.,2.],rtol=0,atol=0)

    def test_separated_intervals_certify_at_least_the_gap(self):
        e=np.array([0.,1.,2.,3.])
        r=distances(e,e,dict(base=[1.,0.,0.]),dict(base=[0.,0.,1.]),'fixed')
        assert_allclose(r['Wasserstein1_q_lower_upper'],[1.,3.],rtol=0,atol=0)
        self.assertEqual(r['maximum_guaranteed_CDF_separation'],1.)


if __name__=='__main__':
    unittest.main()
