"""Refine frozen volume integration on selected warm parents without changing their state."""
import argparse
import json
from math import fsum
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_refinement_transfer_test import source
from early_pole_endpoint_warm_source_cells import local
from early_pole_endpoint_warm_targeted_source_refinement import PARENTS

Q=(0.9082914754476046,*PARENTS,8.280090530852393)


class SelectedVolume(ScalarSubshellChart):
    def __init__(self,bundle,selected,n):
        super().__init__(bundle,SPECS['fine'],2,4,1e-5)
        self.original=self.ns
        self.repeats=np.ones(self.ns,dtype=int)
        self.repeats[selected]=n
        t,w=roots_legendre(n)
        q,wt=[],[]
        for j,(a,b) in enumerate(zip(self.ds.ed[:-1],self.ds.ed[1:],strict=True)):
            if self.repeats[j]==1:
                q.append(self.ds.p[j]); wt.append(self.ds.w[j])
            else:
                v=a+(b-a)*(t+1)/2
                W=w*v*v
                q.extend(v); wt.extend(W*(self.ds.w[j]/fsum(W)))
        self.ds.p=np.asarray(q); self.ds.w=np.asarray(wt)
        self.ns=self.ds.n=len(q)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.initial=self.lift(self.initial)

    def lift(self,y):
        return np.r_[np.repeat(y[:self.original],self.repeats),y[self.original:]]


def run(state,receipt):
    x=state['x_end']
    a,y=unpack(state,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    i=np.searchsorted(ed,Q)
    if list(ed[i])!=list(Q) or min(ed[i])<receipt['warm_floor']:
        raise ArithmeticError('Missing declared diagnostic cells')
    j=i-np.searchsorted(ed,receipt['warm_floor'])
    cfg=load('fine')[0]
    info=source(state,cfg,receipt)
    reference={r['rule']:np.asarray(r['warm_cell_forward_inverse_net'])[j,:]
               for r in receipt['cases'] if r['rule'] in ('Volume32','Volume64')}
    N=info['scale']
    nums=[]
    for n in (32,64,128,256):
        d=SelectedVolume(cfg,i,n)
        if any(not np.array_equal(getattr(d,k),getattr(a,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Charge changed unexpectedly')
        f=d.lift(y)
        pop=np.add.reduceat(d.ds.w*f[:d.ns],np.r_[0,np.cumsum(d.repeats)[:-1]])
        if np.max(abs(pop-a.ds.w*y[:a.ns])/np.maximum(pop,1e-105))>3e-14:
            raise ArithmeticError('Scalar parent populations changed')
        q=d.ds.p
        P=q*q/np.hypot(q,m)
        u,v=[local(d,f,P,x,m,H,ed,k)[i] for k in (4,8)]
        defect=np.max(abs(u-v)/np.maximum(v[:,0,None]+v[:,1,None],1e-100))
        if defect>3e-6:
            raise ArithmeticError('Selected warm radial mesh did not close')
        z=reference.get(f'Volume{n}')
        if z is not None and np.max(abs(v-z)/np.maximum(z[:,0,None]+z[:,1,None],1e-100))>3e-12:
            raise ArithmeticError('Selective lifting failed archived-volume reproducibility')
        nums.append((n,float(defect),v/N))
    high=nums[-1][2][:,2]
    drift=np.max(abs(high-nums[-2][2][:,2]))
    if drift>1e-5:
        raise ArithmeticError(('Higher volume selected-cell meshes underresolved',drift))
    old=info['base'][j]/N
    new=info['split'][j]/N
    improvement=abs(old-high)-abs(new-high)
    return dict(x=x,reference_scale=N,selected_edges=ed[np.stack((i,i+1),axis=1)].tolist(),
                per_order=[dict(order=n,mesh_defect=d,forward_inverse_net=z.tolist()) for n,d,z in nums],
                volume128_to_256_maximum_absolute_change=float(drift),
                full_G4_net=old.tolist(),split_G4_net=new.tolist(),
                improvement_split_minus_full_at_256=improvement.tolist(),
                scope='Only four diagnostic warm cells and two higher finite volume rules. No continuum enclosure or new trajectory.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('receipt',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.state),read(a.receipt)),indent=2,allow_nan=False))
