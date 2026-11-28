"""Test strict matched-sideband admission and preservation of the ordinary path."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import numpy as np
from early_pole_endpoint_cold_native_adjacent_window import FIRST_SIDE_BAND,side


class Window(TestCase):
    def setup_case(self,dir_):
        p=Path(dir_)/'band.jsonl'
        a,b=FIRST_SIDE_BAND
        s=(a+b)/2
        h={'grid':'fine','time_x':.05,'scalar_mass':.01,'omega':.11,'residue':1.,
           'sideband_q_span':[a,b],'sideband_state_sha256':'x'*64,'mesh':[4,6,10,18],
           'mirror_azimuth':True,'total':2}
        z=np.array([[8.,3.,5.],[0.,0.,0.],[.2,.1,.1]])
        z[1]=h['scalar_mass']*z[0]+z[2]
        y=[dict(h,panel=i,bounds=edge,coeff=z.tolist()) for i,edge in enumerate(([a,s],[s,b]))]
        p.write_text(''.join(json.dumps(row)+'\n' for row in y))
        return p,y,dict(h,edges=[a,s,b]),h

    def test_exact_pair_and_bad_state_digest(self):
        with TemporaryDirectory() as t:
            p,y,T,h=self.setup_case(t)
            mod='early_pole_endpoint_cold_native_adjacent_window.topology'
            with patch(mod,side_effect=lambda *_a,**_k:dict(T)):
                e,z=side(p,{},1,h,2)
                self.assertEqual(len(e),3)
                self.assertEqual(z.shape,(2,3,3))
                y[1]['sideband_state_sha256']='wrong'
                p.write_text(''.join(json.dumps(row)+'\n' for row in y))
                with self.assertRaises(ArithmeticError): side(p,{},1,h,2)

    def test_partial_data_fails(self):
        with TemporaryDirectory() as t:
            p,y,T,h=self.setup_case(t)
            p.write_text(json.dumps(y[0])+'\n')
            with patch('early_pole_endpoint_cold_native_adjacent_window.topology',return_value=T):
                with self.assertRaises(ArithmeticError): side(p,{},1,h,2)
