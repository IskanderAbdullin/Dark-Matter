"""Build a genuine selective Gauss2 chart and check its frozen source sum."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_pressure_source_control import cell_weights
from early_pole_endpoint_warm_source_cells import local


class MaskedGauss2(ScalarGauss2):
    def __init__(self,bundle,a,selection):
        super().__init__(bundle,SPECS['fine'],2,4,1e-5,gauss_floor=.01)
        edges=a.ds.ed
        native={(float(q),float(r)):i for i,(q,r) in enumerate(zip(edges[:-1],edges[1:],strict=True))}
        keys=[tuple(q) for q in selection]
        if len(set(keys))!=len(keys) or any(q not in native or q[0]<.01 for q in keys):
            raise ValueError('The mask must consist of distinct complete warm native cells')
        if not np.array_equal(edges,self.ds.ed):
            raise ArithmeticError('The full Gauss and original scalar charts differ in cell edges')
        select=np.zeros(a.ns,dtype=bool)
        select[[native[k] for k in keys]]=True
        # Admissible node pairs already exist in the full warm Gauss chart.
        hi=np.cumsum(self.scalar_repeats)
        lo=np.r_[0,hi[:-1]]
        if np.any(self.scalar_repeats[select]!=2):
            raise ArithmeticError('Selected chart cells have no positive Gauss pair')
        q=[]; w=[]
        for j in range(a.ns):
            if select[j]:
                q.extend(self.ds.p[lo[j]:hi[j]])
                w.extend(self.ds.w[lo[j]:hi[j]])
            else:
                q.append(a.ds.p[j])
                w.append(a.ds.w[j])
        self.scalar_repeats=np.where(select,2,1)
        self.gauss_shells=len(keys)
        self.ds.p=np.array(q)
        self.ds.w=np.array(w)
        self.ds.n=self.ns=len(w)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.initial=self.transfer_step(a.initial)


def check(native,details,mask,n):
    x=native['x_end']
    a,y=unpack(native,False,x)
    candidates=[v for v in mask['rows'] if v['refined_native_cells']==n]
    if details['x']!=x or x not in (mask['train_x'],mask['test_x']) or len(candidates)!=1:
        raise ArithmeticError('Context or mask does not match the source receipt')
    row=candidates[0]
    b=MaskedGauss2(load('fine')[0],a,row['selected_q'])
    z=b.transfer_step(y)
    if any(not np.array_equal(getattr(a,k),getattr(b,k)) for k in ('qd','qp','wd','wp')):
        raise ArithmeticError('The masked chart changed charged quadrature')
    counts=np.add.reduceat(b.ds.w*z[:b.ns],np.r_[0,np.cumsum(b.scalar_repeats)[:-1]])
    if np.max(abs(counts-a.ds.w*y[:a.ns])/np.maximum(a.ds.w*y[:a.ns],1e-105))>2e-14:
        raise ArithmeticError('Selective lift lost native cell count')
    i=b.ds.ed[:-1]>=details['warm_floor']
    m,_,_,H=a.state(x,y)
    P=np.repeat(cell_weights(a.ds.ed,m,details['warm_floor']),b.scalar_repeats)
    select=np.repeat(b.scalar_repeats==2,b.scalar_repeats)
    q=b.ds.p
    P[select]=q[select]**2/np.hypot(q[select],m)
    m1,m2=[local(b,z,P,x,m,H,a.ds.ed,k)[i] for k in (4,8)]
    cc={r['rule']:np.array(r['warm_cell_forward_inverse_net']) for r in details['cases']}
    expected=cc['W1'].copy()
    expected[b.scalar_repeats[i]==2]=cc['Gauss2'][b.scalar_repeats[i]==2]
    div=np.maximum(expected[:,0,None]+expected[:,1,None],1e-100)
    pair=np.max(abs(m1-m2)/div)
    defect=np.max(abs(m2-expected)/div)
    if max(pair,defect)>3e-11:
        raise ArithmeticError(('Real chart did not reproduce the prescribed frozen hybrid',pair,defect))
    ref=cc['Volume64'][:,2]
    value=fsum(m2[:,2]-ref)/fsum(ref)
    target=row['early' if x==mask['train_x'] else 'held_out']['relative_signed_error']
    if abs(target-value)>1e-13:
        raise ArithmeticError('A real chart did not reproduce the held-out pressure contrast')
    return dict(x=x,selected_warm_native_cells=n,scalar_nodes=b.ns,cell_affine_relative_defect=defect,
                cell_radial_order_relative_defect=pair,relative_net_pressure_error=value,
                scope='Actual masked chart for an instantaneous collision source at fixed mass and clock. No masked time evolution or adaptive controller.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('S4',type=Path)
    p.add_argument('details',type=Path)
    p.add_argument('mask',type=Path)
    p.add_argument('--cells',type=int,choices=(3,8,12,22),default=8)
    a=p.parse_args()
    print(json.dumps(check(read(a.S4),read(a.details),read(a.mask),a.cells),indent=2,allow_nan=False))
