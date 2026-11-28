import tempfile
import unittest
import json
from pathlib import Path
import numpy as np
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels
import early_pole_endpoint_cold_targeted_patch_join as join



class RefinedPatch(unittest.TestCase):
    def test_targets_exact_old_panel(self):
        z,c=refine_panels(np.array([1.,2.,3.]),[[2.,2.5,3.]])
        np.testing.assert_array_equal(z,[1.,2.,2.5,3.])
        self.assertEqual(c,{1,2})
        for bad in ([[1.,1.25,3.]],[[2.,2.,3.]],[[2.,2.5,3.],[2.,2.25,3.]]):
            with self.assertRaises(ArithmeticError):
                refine_panels(np.array([1.,2.,3.]),bad)
    
    def test_closure_and_patch_loss_gate(self):
        with tempfile.TemporaryDirectory() as d:
            base,new,out=[Path(d)/x for x in ('base.jsonl','new.jsonl','out.jsonl')]
            edges=[1e-6,2e-6,3e-6,1.5e-5]
            meta=dict(grid='fine', time_x=.05, scalar_mass=1., omega=2., residue=3.,mesh=[2,2,2,2])
            def row(i, n, bounds, v):
                return dict(meta,panel=i,total=n,bounds=bounds,coeff=(v*np.array([[3,1,2],[6,2,4],[3,1,2]])).tolist())
            def dump(path, data):
                path.write_text(''.join(json.dumps(x)+'\n' for x in data))
            dump(base,[row(0,2,edges[:2],2.),row(1,2,[edges[1],edges[-1]],2.)])
            patch=[row(1,3,edges[1:3],1.),row(2,3,edges[2:],1.)]
            plan=dict(meta,total=3,edges=edges,recalculated_panel_indices=[1,2])
            dump(new,patch[:1])
            with self.assertRaises(ArithmeticError):
                join.merge(base,new,plan,out)
            dump(new,patch)
            self.assertEqual(join.merge(base,new,plan,out)['result_panels'],3)
            self.assertEqual(join.reduce(out)[1],[[1e-6,2e-6],[2e-6,3e-6],[3e-6,1.5e-5]])


if __name__=='__main__':
    unittest.main()
