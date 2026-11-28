import unittest
import numpy as np

from early_pole_endpoint_cold_shared_aggregate_separation import compare, support


class SharedAggregateTests(unittest.TestCase):
    def test_support_preserves_correlation_and_spontaneous_offset(self):
        def chart(pairs,offset):
            z=np.array(pairs)
            O=np.array(offset).reshape(3,1)
            return O,np.ones_like(O),[z],[z]
        a=chart([[0.,1.,.5],[1.,0.,.5]],[1.,0.,0.])
        b=chart([[0.,0.,.5]],[1.,0.,0.])
        self.assertEqual(support(a,[1.,1.,0.],0),(2.,2.))
        r=compare(a,b,0)
        self.assertEqual(r['forward_number']['disjoint_gap'],0.)
        self.assertEqual(r['inverse_number']['disjoint_gap'],0.)
        self.assertGreater(r['coupled_three']['guaranteed_linf_fraction_at_least'],.4)


if __name__=='__main__':
    unittest.main()
