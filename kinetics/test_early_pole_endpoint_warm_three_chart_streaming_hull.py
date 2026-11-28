import unittest

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_radiation_clock_streaming_probe import average_speed
from early_pole_endpoint_warm_three_chart_streaming_hull import curvature,hull,kernel


class StreamingEnclosure(unittest.TestCase):
    def test_curvature_against_explicit_derivatives(self):
        t,w=roots_legendre(128)
        u,v=(t+1)/2,w/2
        for a,b,m,phase in ((.01,.016,.00087,990.),(.62,.92,.00087,174.)):
            E=np.hypot(a,m)
            h=np.hypot(b,m)-E
            xi=np.linspace(0.,1.,21)
            q=np.sqrt(a*a+h*xi*(2*E+h*xi))
            q1=h*(E+h*xi)/q
            q2=-h*h*m*m/q**3
            for A in (1.,2.,10.,100.,1000.):
                mass=m*(1+(A-1)*u)
                p=np.hypot(q[:,None],mass)
                s=average_speed(q,m,A)
                v1=(mass**2/p**3)@v
                v2=(-3*mass**2*q[:,None]/p**5)@v
                arg=phase*s[:,None]*u
                F1=phase*np.sin(arg)@(u*v)
                F2=phase**2*np.cos(arg)@(u*u*v)
                D=F2*(v1*q1)**2+F1*(v2*q1*q1+v1*q2)
                self.assertLessEqual(float(np.max(abs(D))),curvature(a,b,m,A,phase)*(1+2e-12))

    def test_two_exactly_feasible_measures_are_enclosed(self):
        a,b,m,A,p=.7,.9,.002,10.,15.
        mu=.38
        E=np.hypot(a,m)
        u=np.hypot(b,m)-E
        outer,near=hull(a,b,m,A,p,mu,96)
        x=np.sqrt(a*a+mu*u*(2*E+mu*u))
        witnesses=[kernel(np.array([x]),m,A,p)[0],kernel(np.array([a,b]),m,A,p)@np.array([1-mu,mu])]
        self.assertLessEqual(outer[0],min(witnesses))
        self.assertGreaterEqual(outer[1],max(witnesses))
        self.assertLessEqual(near[0],witnesses[1])
        self.assertGreaterEqual(near[1],witnesses[1])


if __name__=='__main__':
    unittest.main()
