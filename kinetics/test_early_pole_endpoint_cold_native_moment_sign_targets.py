"""Bound the sign-target planner against enumeration on a tiny nonnested grid."""
from itertools import combinations
import unittest

import numpy as np

from early_pole_endpoint_cold_native_moment_sign_targets import targets
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels


class Targets(unittest.TestCase):
    def test_exact_finite_objective(self):
        e=np.array([1.,2.,4.])
        f=np.array([1.,2.5,3.,4.])
        a=np.array([1.,3.])[:,None,None]*np.ones((1,3,2))
        b=np.array([2.,3.,1.])[:,None,None]*np.ones((1,3,2))
        _,v,w,mask=targets((e,a),(f,b),1.,2)
        self.assertEqual(sum(mask),2)
        scale=np.abs(v)
        self.assertTrue(np.all(scale>0.))
        vals=[]
        for chosen in combinations(range(5),2):
            vals.append(self.upper(e,a[:,0,0],f,b[:,0,0],chosen)/scale[:3])
        self.assertAlmostEqual(max(w/scale),min(max(x) for x in vals),places=12)

    @staticmethod
    def upper(e,a,f,b,chosen):
        rows=[]
        i=0
        for side,(q,z) in enumerate(((e,a),(f,b))):
            let=lambda x: np.array([x,x*x,x/np.hypot(x,1.)])
            total=0.
            for x,y,m in zip(q[:-1],q[1:],z,strict=True):
                left,right=let(x),let(y)
                end=right if side else left
                if i in chosen:
                    end=.5*(end+let((x+y)/2))
                total+=m*end
                i+=1
            rows.append(total/sum(z))
        return rows[1]-rows[0]

    def test_next_split_uses_first_children_as_existing_parents(self):
        q,_=refine_panels(np.array([1.,2.,3.]),[[1.,1.5,2.]])
        v,j=refine_panels(q,[[1.,1.25,1.5]])
        np.testing.assert_array_equal(v,[1.,1.25,1.5,2.,3.])
        self.assertEqual(j,{0,1})
