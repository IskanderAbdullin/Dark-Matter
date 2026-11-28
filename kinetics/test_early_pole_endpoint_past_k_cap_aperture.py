"""Independent closed-form fixtures for full, split and empty parent-cap apertures."""
import unittest

import numpy as np

from early_pole_endpoint_past_k_cap_aperture import aperture


class ParentCapAperture(unittest.TestCase):
    def test_prescribed_cap_cosine(self):
        k,q,scalar,parent,daughter,E=.0015,1.5e-5,.00087,.05,.015,32.
        momentum=np.sqrt((E-parent)*(E+parent))
        for c,expected in ((-1.5,0.),(0.,.5),(1.5,1.)):
            v=np.sqrt((k-q)**2+2*k*q*(1-c))
            z=2*v*momentum-(parent*parent-daughter*daughter)+v*v
            u=z/(np.sqrt(E*E+z)+E)
            w=u+np.hypot(scalar,q)
            a=aperture(k,q,w,scalar,parent,daughter,E)
            self.assertAlmostEqual(a['cosine_upper_before_clipping'],c,places=10)
            self.assertAlmostEqual(a['allowed_cosine_measure_fraction'],expected,places=10)

    def test_invalid_transfer_not_scanned(self):
        with self.assertRaises(ArithmeticError):
            aperture(.0015,1.5e-5,.0005,.00087,.05,.015,32.)


if __name__=='__main__':
    unittest.main()
