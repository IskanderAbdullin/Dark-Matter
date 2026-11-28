"""Analytic and bounded checks independent of the archived trajectories."""
from types import SimpleNamespace
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_warm_observable_quantiles import CUTS, PROBS, observe, parts


class WarmObservableQuantiles(unittest.TestCase):
    def test_single_massless_constant_occupation(self):
        a, b, f = 1e-7, 5., .7
        d = SimpleNamespace(ds=SimpleNamespace(ed=np.array([a,b]),
                                               w=np.array([(b**3-a**3)/(6*np.pi**2)])),ns=1)
        total, q, tail = observe(d,np.array([f]),0.,24)
        powers = np.array([3.,4.,4.,3.])
        expected = (a**powers[:,None]+np.array(PROBS)*(b**powers[:,None]-a**powers[:,None]))**(1/powers[:,None])
        assert_allclose(q,expected,rtol=4e-14,atol=0)
        assert_allclose(total,f*(b**powers-a**powers)/(2*np.pi**2*powers),rtol=4e-14,atol=0)
        z=np.minimum(np.maximum(CUTS,a),b)[:,None]
        assert_allclose(tail,(b**powers-z**powers)/(b**powers-a**powers),rtol=4e-14,atol=0)

    def test_mass_change_cannot_affect_rest_count(self):
        q=np.array([1e-9,2e-8,.2,1.,8.])
        f=np.array([1e5,1e-4,0.,.1])
        a,b=[parts(q[:-1],q[1:],f,m,48) for m in (.001,.1)]
        assert_allclose(a[:,0],b[:,0],rtol=0,atol=0)
        self.assertTrue(np.all(a[f>0,1:]>b[f>0,1:]))
        self.assertTrue(np.all(b[f==0]==0))


if __name__=='__main__':
    unittest.main()
