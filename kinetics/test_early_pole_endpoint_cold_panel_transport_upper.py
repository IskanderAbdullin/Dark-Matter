import unittest
from itertools import product

import numpy as np

from early_pole_endpoint_cold_panel_transport_upper import most


def discrete(loc,wei,other,mass):
    edges=np.union1d(loc,other)
    mid=(edges[:-1]+edges[1:])/2.
    a=np.array([sum(v for u,v in zip(loc,wei) if u<t) for t in mid])
    b=np.array([sum(v for u,v in zip(other,mass) if u<t) for t in mid])
    return np.diff(edges)@abs(a-b)


class PanelWorstCase(unittest.TestCase):
    def test_overlap_forest_matches_exhaustive_endpoint_enumeration(self):
        a=np.array([0.,.2,.64,1.])
        b=np.array([0.,.37,.55,1.])
        p=[.2,.65,.15]
        q=[.41,.18,.41]
        truth=max(discrete(a[np.arange(3)+s[:3]],p,b[np.arange(3)+s[3:]],q)
                  for s in product((0,1),repeat=6))
        r=most(a,b,p,q)
        self.assertAlmostEqual(r['sharp_upper_q'],truth,places=14)
        self.assertAlmostEqual(r['endpoint_witness_q'],truth,places=14)

    def test_interior_split_does_not_exceed_enclosed_extremum(self):
        a=np.array([0.,.6,1.])
        b=np.array([0.,.25,1.])
        p=[.4,.6]
        q=[.66,.34]
        hi=most(a,b,p,q)['sharp_upper_q']
        for s in (.1,.35,.7,.9):
            w=[p[0]*s,p[0]*(1-s),p[1]*s,p[1]*(1-s)]
            v=[q[0]*s,q[0]*(1-s),q[1]*s,q[1]*(1-s)]
            x=[.05,.56,.67,.98]
            y=[.07,.24,.41,.91]
            self.assertLessEqual(discrete(x,w,y,v),hi+3e-15)


if __name__=='__main__':
    unittest.main()
