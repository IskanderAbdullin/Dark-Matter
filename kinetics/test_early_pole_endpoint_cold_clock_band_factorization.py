import unittest
from types import SimpleNamespace
import numpy as np

from early_pole_endpoint_cold_clock_band_factorization import crossed


class ColdClockBandsTests(unittest.TestCase):
    def test_separable_uniform_population(self):
        e=np.array([1.,2.,3.])
        c=SimpleNamespace(ds=SimpleNamespace(ed=e))
        f=[np.ones(2),np.full(2,2.)]
        z=[]
        for amp,n in zip((2.,3.),f,strict=True):
            M=np.zeros((2,3,3))
            M[:,0,0]=amp*(1+n)
            M[:,0,1]=amp*n
            M[:,2,0]=amp*(1+n)
            z.append((e,M))
        A,r=crossed(z,list(zip((c,c),f,strict=True)))
        np.testing.assert_allclose(A[0,1],[12.,8.,12.])
        np.testing.assert_allclose(r['all_bare_kernel_log'],np.full(3,np.log(1.5)))
        np.testing.assert_allclose(r['explicit_population_log'],np.log([1.5,2.,1.5]))
        np.testing.assert_allclose(r['mixed_kernel_population_log_interaction'],0.,atol=1e-14)


if __name__=='__main__':
    unittest.main()
