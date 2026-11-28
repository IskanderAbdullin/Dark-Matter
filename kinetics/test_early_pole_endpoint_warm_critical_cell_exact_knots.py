"""Check the event topology against independently constructed equalities."""
import unittest
from types import SimpleNamespace
import numpy as np
from early_pole_endpoint_warm_critical_cell_exact_knots import knots

class ExactKnots(unittest.TestCase):
    def test_parent_daughter_boundary_collision_is_present(self):
        x,m,ratio=.1,.002,.5
        p=np.linspace(.1,5.,44)
        f=np.linspace(.3,4.,31)
        d=SimpleNamespace(p_edges=p,d_edges=f,ratio=ratio)
        e=np.hypot(p[:,None],x)-np.hypot(f[None,:],x*ratio)
        j=e[(e>.4)&(e<1.2)]
        q=np.sqrt((j-m)*(j+m))
        k=knots(d,x,m,.4,1.2)
        self.assertEqual(k[0],.4)
        self.assertEqual(k[-1],1.2)
        want=q[(q>.4)&(q<1.2)]
        self.assertEqual(np.count_nonzero(np.isin(want,k)),len(want))
        self.assertTrue(np.all(np.diff(k)>0.))
