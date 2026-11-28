"""Tamper checks for the warm nodal receipt without rerunning its evolution."""
import copy
from types import SimpleNamespace
import unittest

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import S,X,verify


def example():
    d=SimpleNamespace(ns=3,nf=1,np=1,n=5,old_ns=2,scalar_repeats=np.array([1,2]),
                      ds=SimpleNamespace(w=np.array([1.,.1,.2])),
                      state=lambda x,y:(2e-5,None))
    chart=dict(grid='fine',charge_div=2,scalar_div=4,cut=1e-5,gauss_floor=.01,reference_x=X,
               nodes=[3,1,1],scalar_cell_count=2,scalar_nodal_rule='gauss-volume-v1')
    target=.0075
    steps=2048
    y=[0.1,0.2,0.4,0.2,0.6]
    v=np.multiply(d.ds.w,y[:d.ns])/S
    cells=[v[0],v[1]+v[2]]
    ends=X*np.exp(np.log(target/X)*np.arange(1,9)/8)
    rows=[dict(x=t,scalar_yield=sum(v),pair_number_error=0.,continuity_residual=0.,
               max_charged=.6,max_limited_extent=0.) for t in ends]
    state=dict(chart,x_end=target,steps=2*steps,nodal_occupation=y)
    report=dict(chart,starting_x=X,x_end=target,coarse_steps=steps,
                checkpoint_trajectories=[rows,copy.deepcopy(rows)],scalar_cell_yields=cells,
                scalar_pole_mass=2e-5,relative_time_drift=0.)
    return report,state,d,X,target,steps


class GaussWarmReceipt(unittest.TestCase):
    def test_accepts_consistent_data(self):
        self.assertGreater(verify(*example())['scalar_yield'],0)

    def test_wrong_cell_counts_despite_same_total(self):
        r,*rest=example()
        r['scalar_cell_yields'][0]*=1.01
        r['scalar_cell_yields'][1]*=.99
        with self.assertRaises(ArithmeticError):
            verify(r,*rest)

    def test_nonterminal_mislabeled_clock(self):
        r,*rest=example()
        r['checkpoint_trajectories'][0][0]['x']*=1.0000001
        with self.assertRaises(ArithmeticError):
            verify(r,*rest)

    def test_finite_reported_drift(self):
        r,*rest=example()
        r['relative_time_drift']=float('nan')
        with self.assertRaises(ArithmeticError):
            verify(r,*rest)


if __name__=='__main__':
    unittest.main()
