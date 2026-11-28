import unittest

from finite_mass_qed_native_past_cutoff_check import quantiles
import numpy as np


class NativeCutoff(unittest.TestCase):
    def test_quantile_brackets_are_native_not_interpolated(self):
        q=np.array([0.,4.,5.,6.])
        p=np.array([0.,.7,.9,1.])
        self.assertEqual(quantiles(q,p,[.5,.9,.99]),[[0.,4.],[4.,5.],[5.,6.]])

    def test_negative_prefix_step_is_inadmissible(self):
        with self.assertRaises(ArithmeticError):
            quantiles(np.arange(4),[0.,.6,.5,1.],[.5])


if __name__=='__main__':
    unittest.main()
