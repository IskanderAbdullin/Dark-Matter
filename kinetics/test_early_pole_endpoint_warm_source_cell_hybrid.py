"""Check that held-out info does not enter the frozen mask."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from early_pole_endpoint_warm_source_cell_hybrid import run


def card(x,errors):
    ref=np.ones(len(errors))
    z=np.zeros_like(ref)
    return dict(x=x,warm_native_edges=list(range(len(ref)+1)),ref_net_pressure=float(len(ref)),
                mesh_reference=dict(global_L1_over_total=True),
                cases=[dict(rule=k,warm_cell_forward_inverse_net=np.stack((n,z,n),axis=-1).tolist())
                       for k,n in (('W1',ref+errors),('Gauss2',ref),('Volume64',ref))])


class FrozenMask(unittest.TestCase):
    def test_strictly_early_ranking(self):
        a=np.zeros(22)
        a[[4,9,3]]=[.5,.3,.1]
        b=np.zeros(22)
        b[[2,4,9,3]]=[2.,.5,.3,.1]
        out=run(card(.01,a),card(.05,b))
        r=next(r for r in out['rows'] if r['refined_native_cells']==3)
        self.assertEqual(r['selected_q'],[[4,5],[9,10],[3,4]])
        assert_allclose(r['early']['relative_L1_residual'],0.,atol=1e-16)
        assert_allclose(r['held_out']['relative_signed_error'],2./22,rtol=1e-15)

    def test_requires_same_edges(self):
        a=card(.01,np.zeros(22))
        b=card(.05,np.zeros(22))
        b['warm_native_edges'][-1]+=1
        with self.assertRaises(ArithmeticError):
            run(a,b)


if __name__=='__main__':
    unittest.main()
