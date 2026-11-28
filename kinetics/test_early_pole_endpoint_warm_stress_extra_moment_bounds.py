import unittest
import numpy as np

from early_pole_endpoint_warm_stress_extra_moment_bounds import local
from early_pole_endpoint_warm_stress_conditional_hull import values


class ExtraMomentBoundsTests(unittest.TestCase):
    def test_local_tail_bounds(self):
        q=np.array([1.001,1.1,1.3,1.998])
        w=np.array([.1,.5,.1,.3])
        for m in (1.,3.,10.,100.):
            y=w@values(q*q,m)
            for i,v in enumerate(y,1):
                for k in range(7):
                    center,err,lo,hi=local(1.,2.,q,w,m,i,k)
                    self.assertLessEqual(max(center-err,lo),v+1e-15)
                    self.assertGreaterEqual(min(center+err,hi),v-1e-15)


if __name__=='__main__':
    unittest.main()
