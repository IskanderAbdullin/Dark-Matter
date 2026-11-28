"""Compare interrupted and uninterrupted actual short endpoint runs."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive,write_drive_checkpoint,drive_checkpoint_digest
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_stencil_remap_pilot import X,load


class Pause(Exception):
    pass


class Restartable(unittest.TestCase):
    def test_exact_short_actual_restart(self):
        d=WarmGauss(load('fine')[0],8)
        n=16
        T=X*1.002
        opt=dict(stepper=step)
        a,y=drive(d,T,n,**opt)
        prov='a'*64
        with TemporaryDirectory() as k:
            f=Path(k)/'state.json'
            def preserve(p):
                write_drive_checkpoint(f,p)
                if p['completed']==6:
                    raise Pause()
            with self.assertRaises(Pause):
                drive(d,T,n,checkpoint=preserve,checkpoint_stride=3,provenance=prov,**opt)
            cp=json.loads(f.read_text())
            b,z=drive(d,T,n,restart=cp,provenance=prov,**opt)
            self.assertEqual(a,b)
            np.testing.assert_array_equal(y,z)
            with self.assertRaises(ArithmeticError):
                drive(d,T*1.00001,n,restart=cp,provenance=prov,**opt)
            with self.assertRaisesRegex(ArithmeticError,'Corrupt'):
                drive(d,T,n,restart=dict(cp,occupation=[0.]*d.n),provenance=prov,**opt)
            bad=dict(cp,trace=cp['trace'][:-1])
            bad['payload_sha256']=drive_checkpoint_digest(bad)
            with self.assertRaisesRegex(ArithmeticError,'unresolved'):
                drive(d,T,n,restart=bad,provenance=prov,**opt)

    def test_completed_arm_requires_no_new_collision_calls(self):
        from unittest.mock import patch
        d=WarmGauss(load('fine')[0],8)
        T=X*1.001
        keep=[]
        a,y=drive(d,T,8,stepper=step,checkpoint=keep.append,checkpoint_stride=8,provenance='b'*64)
        self.assertEqual((keep[-1]['completed'],keep[-1]['target_steps']),(8,8))
        with patch.object(d,'chart',side_effect=AssertionError('completed arm was re-evolved')):
            b,z=drive(d,T,8,stepper=step,restart=keep[-1],provenance='b'*64)
        self.assertEqual(a,b)
        np.testing.assert_array_equal(y,z)
