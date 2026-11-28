import unittest

import numpy as np

from early_pole_endpoint_cold_profile_interval_transport import fixed,repair_flow,transport


class IntervalTransport(unittest.TestCase):
    def test_panels_touching_at_an_edge_cost_zero_in_support_relaxation(self):
        a,b=[fixed({'base':v}) for v in ([.8,.2],[.2,.8])]
        p=transport(a,b,np.array([0.,.5,1.]),np.array([0.,.5,1.]),'fixed',independent_bracket=True)
        lo,hi=p['independent_cost_lower_upper_q']
        self.assertEqual(lo,0.)
        self.assertLess(hi,1e-12)

    def test_disjoint_nonadjacent_panels_need_positive_transport(self):
        a,b=[fixed({'base':v}) for v in ([1.,0.,0.],[0.,0.,1.])]
        z=np.array([0.,.2,.7,1.])
        p=transport(a,b,z,z,'fixed',independent_bracket=True)
        lo,hi=p['independent_cost_lower_upper_q']
        self.assertLess(lo,.5)
        self.assertGreater(hi,.5)
        self.assertLess(hi-lo,1e-12)

    def test_projection_respects_both_corrected_marginals(self):
        r,e=repair_flow(np.array([[.15,-.1],[.7,.35]]),np.array([.4,.6]),np.array([.7,.3]))
        self.assertLess(e,2e-14)
        np.testing.assert_allclose(r.sum(axis=0),[.7,.3])
        np.testing.assert_allclose(r.sum(axis=1),[.4,.6])
        self.assertGreaterEqual(r.min(),0.)


if __name__=='__main__':
    unittest.main()
