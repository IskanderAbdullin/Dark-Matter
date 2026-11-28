import unittest

import numpy as np
from scipy.linalg import expm

from fd_portal_direct_comoving_variational_jacobian import diff,full
from fd_portal_direct_comoving_variational_propagator import ordered


class OrderedComoving(unittest.TestCase):
    def test_collisionless_radiation_is_only_baseline_dilution(self):
        class Radiation:
            lo,hi=.05,2.
            eos=None
            def __init__(self):
                self.eos=self
            def state(self,T):
                rho=2.*T**4
                return rho,rho/3.,(8./3.)*T**3,8.*T**3
            def bath(self,u):
                eta,l=u
                T=np.exp(l)
                n=np.exp(eta)*T**3
                return dict(n=n,rho=3.*n*T,pressure=n*T,dn_deta=n,chemical_energy=3.*T,cv_fixed_number=3.*n)
            def collision(self,*unused):
                return np.zeros(2),0.
        f=full(Radiation(),[.6,.04,np.log(.57)])
        np.testing.assert_allclose(f(np.array([1.24,.8,.91])),[[0.,-.8,-.91],[0.,0.,0.]],rtol=2e-13,atol=1e-14)
        end=.24
        y=lambda x:np.exp(-x*np.array([0.,1.,1.]))
        R,_,_=ordered(lambda x,z:f(z),y,np.linspace(0.,end,9),3e-4)
        np.testing.assert_allclose(R,np.diag(y(end)),rtol=2e-11,atol=2e-13)
        np.testing.assert_allclose(R/y(end)[:,None],np.eye(3),rtol=2e-11,atol=2e-13)
        # A constant physical collision vector changes only through H when U varies.
        class FixedSources(Radiation):
            def collision(self,*unused):
                return np.array([2e-21,-3e-21]),0.
        g=full(FixedSources(),[.6,.04,np.log(.57)])
        C=g(np.ones(3))[1]
        A=diff(g,3e-4)[1]
        np.testing.assert_allclose(A[:,:2],0.,atol=2e-14)
        np.testing.assert_allclose(A[:,2],-.5*C,rtol=2e-10,atol=2e-14)
        np.testing.assert_array_equal(A[2],0.)

    def test_noncommuting_with_nonstationary_normalization(self):
        a,c,d,b,theta=.6, .3, .7, 1.1, .24
        def matrix(x):
            return np.array([[-a,b*x,0.],[0.,-c,0.],[0.,0.,-d]])
        def rates(x,v):
            M=matrix(x)
            return np.array([np.diag(M.diagonal())@v,(M-np.diag(M.diagonal()))@v])
        def center(x):
            return np.array([1.+.11*x,np.exp(-.8*x),.98-.2*x])
        P,error,trace=ordered(rates,center,np.linspace(0.,theta,13),2e-4)
        mean=a-c
        target=np.diag(np.exp(-theta*np.array([a,c,d])))
        target[0,1]=np.exp(-a*theta)*b*((mean*theta-1.)*np.exp(mean*theta)+1.)/mean**2
        np.testing.assert_allclose(P,target,rtol=3e-10,atol=2e-13)
        self.assertLess(error,2e-12)
        self.assertAlmostEqual(trace,-(a+c+d)*theta,places=12)
        averaged=expm(theta*matrix(theta/2.))
        self.assertGreater(abs(P[0,1]-averaged[0,1]),1e-4)


if __name__=='__main__':
    unittest.main()
