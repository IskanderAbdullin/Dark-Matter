from itertools import product
import unittest

from early_pole_endpoint_cold_omitted_response_gate import break_order,tv_interval,worst_rank_gap


class OmittedResponse(unittest.TestCase):
    def test_exact_exterior_limits_for_unequal_fractions(self):
        lo,hi=tv_interval([.6,.4],[.4,.6],.1,.2)
        self.assertAlmostEqual(lo,.22)
        self.assertAlmostEqual(hi,.32)

    def test_worst_case_common_fraction_has_a_sharp_threshold(self):
        a,b,c=[1.,0.],[.7,.3],[.6,.4]
        e=break_order(.3,.1)
        self.assertAlmostEqual(e,1./6.)
        for delta,s in ((-.05,-1),(0.,0),(.05,1)):
            t=e+delta
            u=tv_interval(a,b,t,t)[0]
            v=tv_interval(b,c,t,t)[1]
            self.assertAlmostEqual(v-u,1.2*delta)
            self.assertEqual((v-u>1e-15)-(v-u<-1e-15),s)

    def test_adversarial_unequal_caps_are_bounded_by_shared_witness(self):
        a,b,c=[.55,.15,.30],[.12,.53,.35],[.18,.43,.39]
        for cap in (.01,.15,.45,.94):
            worst=worst_rank_gap(a,b,c,cap)
            for e,f,g in product((0.,cap/5.,cap/2.,cap),repeat=3):
                x=tv_interval(a,b,e,f)[0]
                y=tv_interval(b,c,f,g)[1]
                self.assertLessEqual(y-x,worst+1e-14)


    def test_invalid_mass_fraction_rejected(self):
        with self.assertRaises(ArithmeticError):
            tv_interval([1.],[1.],-.1,.2)


if __name__=='__main__':
    unittest.main()
