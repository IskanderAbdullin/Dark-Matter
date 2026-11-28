"""Explicit bin CDF construction and the joint release bound."""
import unittest
from itertools import product
import numpy as np
from early_pole_endpoint_cold_spatial_shape_release import A, B, adversary

class ShapeRelease(unittest.TestCase):
    def test_joint_release_bound_on_cube_vertices(self):
        p=np.array([.22,.12,.15,.17])
        for x in product((0.,1.), repeat=5):
            d=p*(np.array(x[:4])-x[-1])
            self.assertLessEqual(abs(A@d)+abs(B@d),max(p)+1e-16)

    def test_kernel_can_dominate_even_when_population_previously_does(self):
        e=np.array([0.,1.,2.,3.])
        p=np.array([[.2,.3,.5],[.23,.3,.47],[.2,.3,.5],[.23,.3,.47]])
        f=adversary(e,p)
        a,b=f(0.)
        self.assertGreater(a,b)
        a,b=f(.9)
        self.assertLess(a,b)
