import unittest
import numpy as np

from early_pole_endpoint_cold_panel_reweight_envelope import panel_excess
from early_pole_endpoint_cold_shared_reference_scan import feasibility, max_feasible_common_partition


class SharedReferenceTests(unittest.TestCase):
    def test_maximal_merge_recovers_both_reference_moments(self):
        coarse=np.array([1.,2.,3.,4.])
        fine=np.array([1.,1.5,2.,3.,4.])
        mass=1.
        base=panel_excess(coarse[:-1],coarse[1:],mass,48)
        shifted=base.copy()
        shifted[0]=panel_excess(np.array([1.5]),np.array([2.]),mass,48)[0]
        n=np.ones(3)
        self.assertEqual(len(feasibility(coarse,coarse,mass,shifted,n)),1)
        self.assertEqual(len(feasibility(fine,coarse,mass,shifted,n)),0)
        i=max_feasible_common_partition(coarse,[(n,base),(n,shifted)],[coarse,fine],mass)
        np.testing.assert_array_equal(i,[0,2,3])
        for means in (base,shifted,0.4*base+0.6*shifted):
            value=np.add.reduceat(means,i[:-1])/np.diff(i)
            self.assertEqual(feasibility(coarse,coarse[i],mass,value,np.diff(i)),[])
            self.assertEqual(feasibility(fine,coarse[i],mass,value,np.diff(i)),[])


if __name__=='__main__':
    unittest.main()
