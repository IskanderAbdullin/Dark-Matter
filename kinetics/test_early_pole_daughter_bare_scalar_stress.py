"""Analytic scalar stress checks without collision or archive inputs."""
import unittest
import numpy as np
from early_pole_daughter_bare_scalar_stress_control import state


class ScalarTrace(unittest.TestCase):
    def test_zero_momentum(self):
        p=np.array([0.])
        n=np.array([1.])
        # m=2, bare=1.5; 2(4-2.25) = lambda/2.
        r=state(p,n,7.,1.5)
        self.assertAlmostEqual(r['effective_gap'],2.,places=12)
        self.assertAlmostEqual(r['energy']-3*r['pressure'],1.5**2/2,places=12)
        c=state(p,n,7.,0.)
        self.assertAlmostEqual(c['effective_gap'],3.5**(1/3),places=12)
        self.assertAlmostEqual(c['pressure_to_energy'],1/3,places=12)


if __name__=='__main__':
    unittest.main()
