"""Independent drift checks with a controlled analytic reservoir."""
import unittest

import numpy as np

from fd_portal_tabulated_entropy_split_expansion import Stage, hubble


class Ideal:
    lo,hi=np.log((.009,1.05))
    def state(self,t):
        c=10.75*np.pi**2/30
        return c*t**4,c*t**4/3,4*c*t**3/3,4*c*t**3


class ExpansionInvariants(unittest.TestCase):
    def test_adiabatic_pressure_work_and_count(self):
        a=np.array([.05,.3,np.log(.058)])
        g=Stage(Ideal(),None,112)
        initial=g.bath(a[1:])
        h=hubble(g.eos.state(a[0])[0],initial['rho'])
        x,b,w,d=g.drift(0.,a,.005/h)
        final=g.bath(b[1:])
        self.assertGreater(w,0.)
        self.assertLess(d,7e-13)
        self.assertAlmostEqual(g.energy(x,b)/g.energy(0.,a)+w/g.energy(0.,a),1.,places=12)
        self.assertAlmostEqual(g.entropy(x,b)/g.entropy(0.,a),1.,places=12)
        self.assertAlmostEqual(np.exp(3*x)*final['n']/initial['n'],1.,places=12)

    def test_direct_drift_cooling_given_scale_change(self):
        x=.008
        a=np.array([.3,-.5,np.log(.35)])
        g=Stage(Ideal(),None,112)
        b=g.drift_chart(a,x)
        self.assertAlmostEqual(b[0]/a[0],np.exp(-x),places=13)
        self.assertLess(b[2],a[2])
        self.assertAlmostEqual(g.entropy(x,b)/g.entropy(0.,a),1.,places=12)
