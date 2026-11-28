import unittest

from early_pole_endpoint_cold_native_aperture_rectangle import threshold


class ParentAperture(unittest.TestCase):
    def test_independent_threshold_derivatives(self):
        m1,m2,t,v=.05,.003,-.00143,.0015
        e,dt,dv=threshold(m1,m2,t,v)
        self.assertGreater(e,15.)
        h=1e-9
        numeric_t=(threshold(m1,m2,t+h,v)[0]-threshold(m1,m2,t-h,v)[0])/(2*h)
        numeric_v=(threshold(m1,m2,t,v+h)[0]-threshold(m1,m2,t,v-h)[0])/(2*h)
        self.assertAlmostEqual(dt/numeric_t,1.,places=8)
        self.assertAlmostEqual(dv/numeric_v,1.,places=8)

    def test_nonpast_and_lightcone_inputs_fail_closed(self):
        for t,v in ((.001,.002),(-.002,.002)):
            with self.assertRaises(ArithmeticError):
                threshold(.05,.003,t,v)


if __name__=='__main__':
    unittest.main()
