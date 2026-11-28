import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from early_pole_endpoint_gauss_warm_comparison import chart_for
from early_pole_endpoint_targeted_warm_history import build,run


class TargetedWarmHistoryTests(unittest.TestCase):
    def test_metadata_and_positive_lift(self):
        d,s=build()
        kind,c,rep=chart_for(s)
        self.assertEqual(kind,'targeted-2xG4')
        self.assertEqual(s['nodes'],[219,448,488])
        self.assertEqual(np.count_nonzero(rep==8),2)
        np.testing.assert_array_equal(d.ds.w,c.ds.w)
        self.assertTrue(np.all(c.ds.w>0))

    def test_restart_rejects_different_split(self):
        d,s=build()
        s['targeted_rule']['split_fraction']=.45
        with TemporaryDirectory() as z:
            a=Path(z)/'x01.json'
            a.write_text(json.dumps(s))
            with self.assertRaisesRegex(ArithmeticError,'Mismatched targeted warm restart'):
                run(.02,1024,a)


if __name__=='__main__':
    unittest.main()
