"""Compare two-cell order increases with two-cell splitting at fixed occupations."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
from early_pole_endpoint_warm_source_cells import local
from early_pole_endpoint_warm_targeted_source_refinement import PARENTS,Targeted


class SelectiveOrder(WarmGauss):
    def __init__(self,bundle,order):
        if order not in (5,6):
            raise ValueError('Keep base order four and increase only the fixed two cells')
        super().__init__(bundle,4)
        new=WarmGauss(bundle,order)
        ia=np.r_[0,np.cumsum(self.repeats)]
        ib=np.r_[0,np.cumsum(new.repeats)]
        chosen=np.isin(self.ds.ed[:-1],PARENTS)
        if sum(chosen)!=2 or np.any(self.repeats[chosen]!=4) or not np.array_equal(self.ds.ed,new.ds.ed):
            raise ArithmeticError('Selective rule changed the archived parent structure')
        initial=np.r_[self.initial[ia[:-1]],self.initial[ia[-1]:]]
        q,w=[],[]
        for j,k in enumerate(chosen):
            d,i=(new,ib) if k else (self,ia)
            s=slice(i[j],i[j+1])
            q.extend(d.ds.p[s])
            w.extend(d.ds.w[s])
        self.repeats=np.where(chosen,order,self.repeats)
        self.ds.p=np.asarray(q); self.ds.w=np.asarray(w)
        self.ns=self.ds.n=len(q)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.initial=self.transfer_step(initial)


def run(state,reference):
    x=state['x_end']; a,y=unpack(state,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    j=np.searchsorted(ed,reference['warm_floor'])
    e=ed[j:]
    if (x,m,H)!=(reference['x'],reference['mass_held'],reference['clock_held']) or e.tolist()!=reference['warm_native_edges']:
        raise ArithmeticError('Frozen snapshot and reference differ')
    if reference['mesh_reference']['global_L1_fraction']>4e-5:
        raise ArithmeticError('Refinement needs an accepted reference')
    N=a.ds.w*y[:a.ns]
    C={r['rule']:np.asarray(r['warm_cell_forward_inverse_net']) for r in reference['cases']}
    ref=C['Volume64']; scale=fsum(ref[:,2])
    if scale<=0 or abs(scale/reference['ref_net_pressure']-1)>2e-14:
        raise ArithmeticError('Reference normalization changed')
    cfg=load('fine')[0]
    mask=np.isin(e[:-1],PARENTS)
    tabs=[]
    for label,d in [('G4',WarmGauss(cfg,4)),('G5',WarmGauss(cfg,5)),('G6',WarmGauss(cfg,6)),
                    ('selective5',SelectiveOrder(cfg,5)),('selective6',SelectiveOrder(cfg,6)),('split4',Targeted(cfg))]:
        u=d.transfer_step(y)
        ps=np.add.reduceat(d.ds.w*u[:d.ns],np.r_[0,np.cumsum(d.repeats)[:-1]])
        if np.max(abs(ps-N)/np.maximum(N,1e-100))>3e-14 or any(
                not np.array_equal(getattr(a,k),getattr(d,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Changed archived populations or charge chart')
        q=d.ds.p; R=np.where(q>=e[0],q*q/np.hypot(q,m),0.)
        lo,hi=[local(d,u,R,x,m,H,ed,k)[j:] for k in (4,8)]
        dlt=float(np.max(abs(lo-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100)))
        if dlt>3e-6:
            raise ArithmeticError('Warm source not on its paired internal meshes')
        tabs.append((label,d.ns,hi,dlt))
    base=tabs[0][2]
    for k,(_,_,xv,_) in enumerate(tabs[3:],3):
        tgt=base.copy()
        if k<5:
            tgt[mask]=tabs[k-2][2][mask]
            err=np.max(abs(xv-tgt)/np.maximum(abs(tgt[:,0,None])+abs(tgt[:,1,None]),1e-100))
        else:
            err=np.max(abs(xv[~mask]-tgt[~mask])/np.maximum(abs(tgt[~mask,0,None])+abs(tgt[~mask,1,None]),1e-100))
        if err>3e-13:
            raise ArithmeticError('Whole hybrid chart contradicts cellwise replacement')
    out=[]
    for label,n,v,e0 in tabs:
        D=(v[:,2]-ref[:,2])/scale
        W=(v[:,2]-C['Volume32'][:,2])/scale
        out.append(dict(rule=label,nodes=n,total_signed_error=fsum(D),L1_vs_64=fsum(abs(D)),
                        L1_vs_32_in_64_scale=fsum(abs(W)),selected_cell_L1_vs_64=fsum(abs(D[mask])),
                        internal_mesh_defect=e0))
    return dict(x=x,warm_edges=e.tolist(),selected_cell_lower_edges=list(PARENTS),scores=out,
                scope='Four actual frozen reference snapshots with the same preselected cells. No selective-order dynamics or continuum claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('reference',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.state),read(a.reference)),indent=2,allow_nan=False))
