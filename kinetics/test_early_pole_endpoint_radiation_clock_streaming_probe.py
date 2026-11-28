import unittest
import numpy as np
from early_pole_endpoint_radiation_clock_streaming_probe import average_speed,clock_check,ratio_block


class RadiationClockProbeTests(unittest.TestCase):
    def test_free_path_limits_and_independent_clock(self):
        q=np.geomspace(1e-6,30.,32)
        self.assertTrue(np.array_equal(average_speed(q,0.,10.),np.ones_like(q)))
        self.assertTrue(np.array_equal(average_speed(q,.003,1.),q/np.hypot(q,.003)))
        self.assertLess(clock_check(q,.003,1000.),1e-12)
        np.testing.assert_allclose(average_speed(q,1e11,10.)/(q/1e11),np.log(10)/9,rtol=5e-15)

    def test_symmetric_split_cancels_when_measures_identical(self):
        q=np.array([.1,.5])
        w=np.array([.8,.2])
        a,_,_=ratio_block([(q,w)]*3,.3)
        np.testing.assert_array_equal(a[9:],np.zeros(8))


if __name__=='__main__':
    unittest.main()
