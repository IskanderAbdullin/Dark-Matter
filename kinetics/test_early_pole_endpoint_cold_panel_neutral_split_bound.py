"""Check the conservative correlated bound against independently enumerated vertices."""
from itertools import product
import unittest

import numpy as np

from early_pole_endpoint_cold_panel_neutral_split_bound import bounds


class ColdSplitBound(unittest.TestCase):
    def test_vertices_and_interior_are_inside_the_envelope(self):
        offset=np.full(4,.3)
        choices=[np.array([[.33,.30,.39,.35],[.37,.40,.31,.35]]),
                 np.array([[.30,.42,.33,.40],[.40,.28,.37,.30]])]
        a,b,m=1.05e-8,1.e-8,.001
        packed=bounds(offset,choices,a,b,m)
        lo,hi=packed['conditional_relative_split_interval']
        for s,t in product((0.,.25,.5,.75,1.),repeat=2):
            x=offset+s*choices[0][0]+(1-s)*choices[0][1]+t*choices[1][0]+(1-t)*choices[1][1]
            B,D,X,Y=x
            r=(a*X*D-b*B*Y)/(B*(m*D+b*Y))
            self.assertLessEqual(lo,r)
            self.assertGreaterEqual(hi,r)
        self.assertAlmostEqual(packed['baseline_relative_split'],(a-b)/(m+b),places=18)


if __name__=='__main__':
    unittest.main()
