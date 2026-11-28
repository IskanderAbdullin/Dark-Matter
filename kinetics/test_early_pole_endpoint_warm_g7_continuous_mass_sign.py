"""Test the analytic interpolation remainder independently at off-grid energies."""
import unittest

import numpy as np
from early_pole_endpoint_warm_g7_continuous_mass_sign import tighter_error
from early_pole_endpoint_warm_stress_conditional_hull import values


class ContinuousMassSign(unittest.TestCase):
    def test_remainder_bounds(self):
        for a,b in ((.01,.013),(.2,.3),(4.,5.8)):
            for m0 in (.00026,.08):
                for mass in (0.,m0,.005,.1,.3):
                    for n in (4,9):
                        E=np.linspace(np.hypot(a,m0),np.hypot(b,m0),n+1)
                        f=lambda e: values((e-m0)*(e+m0),mass)
                        u=f(E[:-1])
                        v=f(E[1:])
                        for t in (.11,.23,.5,.77,.89):
                            mid=f((1-t)*E[:-1]+t*E[1:])
                            err=abs(mid-((1-t)*u+t*v))
                            guard=80*np.finfo(float).eps*np.maximum(abs(u),abs(v))
                            self.assertTrue(np.all(err<=tighter_error(a,b,mass,m0,n)+guard))

    def test_pointwise_decrease_in_mass_and_increase_in_q(self):
        z=np.geomspace(1e-12,100.,31)
        for m in (0.,1e-4,.1,1.):
            self.assertTrue(np.all(values(z,m+1e-3)<values(z,m)))
            self.assertTrue(np.all(np.diff(values(z,m),axis=0)>0.))
