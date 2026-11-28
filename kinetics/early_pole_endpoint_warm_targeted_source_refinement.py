"""Compare targeted within-cell quadrature with uniform Gauss4 and Gauss5."""
import argparse
import json
from math import fsum
from pathlib import Path
import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss,volume_gauss
from early_pole_endpoint_warm_source_cells import local

PARENTS=(3.9636789434531003,5.728841112042199)


class Targeted(WarmGauss):
    def __init__(self,bundle,parents=PARENTS):
        super().__init__(bundle,4)
        edges=self.ds.ed
        i=np.searchsorted(edges,parents)
        if np.any(i>=len(self.repeats)) or len(set(i))!=len(parents) or list(edges[i])!=list(parents) or np.any(self.repeats[i]!=4):
            raise ArithmeticError('Expected exact predeclared warm parent boundaries')
        target=set(i.tolist())
        start=np.r_[0,np.cumsum(self.repeats)]
        q,w=[],[]
        for j,(a,b) in enumerate(zip(edges[:-1],edges[1:],strict=True)):
            sl=slice(start[j],start[j+1])
            if j not in target:
                q.extend(self.ds.p[sl])
                w.extend(self.ds.w[sl])
                continue
            mid=(a+b)/2
            n,v,_=volume_gauss(np.array([a,mid]),np.array([mid,b]),4)
            q.extend(n.ravel())
            w.extend(v.ravel()*(fsum(self.ds.w[sl])/fsum(v.ravel())))
            self.repeats[j]=8
        self.ds.p=np.asarray(q)
        self.ds.w=np.asarray(w)
        self.ns=self.ds.n=len(q)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        # WarmGauss initial already contains nodal values, so collapse the old
        # constant lift before expanding on the selective test chart again.
        a=np.r_[self.initial[start[:-1]],self.initial[start[-1]:]]
        self.initial=self.transfer_step(a)


def run(state,ref):
    x=state['x_end']
    a,y=unpack(state,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    floor=ref['warm_floor']
    masked=ed[:-1]>=floor
    edw=ed[np.searchsorted(ed,floor):]
    if (x,m,H)!=(ref['x'],ref['mass_held'],ref['clock_held']) or edw.tolist()!=ref['warm_native_edges']:
        raise ArithmeticError('Changed state or warm reference')
    raw={c['rule']:np.asarray(c['warm_cell_forward_inverse_net']) for c in ref['cases']}
    R=raw['Volume64']
    S=fsum(R[:,2])
    if S<=0 or abs(S/ref['ref_net_pressure']-1)>2e-14 or ref['mesh_reference']['global_L1_fraction']>4e-5:
        raise ArithmeticError('Unresolved finite warm reference')
    cfg=load('fine')[0]
    plots=[]
    for g in (WarmGauss(cfg,4),WarmGauss(cfg,5),Targeted(cfg)):
        f=g.transfer_step(y)
        n=np.add.reduceat(g.ds.w*f[:g.ns],np.r_[0,np.cumsum(g.repeats)[:-1]])
        if np.max(abs(n-a.ds.w*y[:a.ns])/np.maximum(n,1e-105))>3e-14:
            raise ArithmeticError('Changed native scalar count')
        if any(not np.array_equal(getattr(g,k),getattr(a,k)) for k in ('wp','wd','qp','qd')):
            raise ArithmeticError('Changed charge chart')
        q=g.ds.p
        P=np.where(q>=floor,q*q/np.hypot(q,m),0.)
        lo,hi=[local(g,f,P,x,m,H,ed,k)[masked] for k in (4,8)]
        err=float(np.max(abs(lo-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100)))
        if err>3e-6:
            raise ArithmeticError('Unresolved paired source')
        plots.append((g.ns,hi,err))
    keep=~np.isin(edw[:-1],PARENTS)
    drift=np.max(abs(plots[2][1][keep]-plots[0][1][keep])/np.maximum(np.max(abs(plots[0][1][keep]),axis=1)[:,None],1e-100))
    if drift>3e-13:
        raise ArithmeticError(('Refinement unexpectedly changed untreated cells',drift))
    names=('uniform_G4','uniform_G5','targeted_2xG4')
    summary=[]
    for tag,(n,j,e) in zip(names,plots,strict=True):
        q=(j[:,2]-R[:,2])/S
        summary.append(dict(rule=tag,scalar_nodes=n,reference_signed_error=fsum(q),
                            reference_l1_error=fsum(abs(q)),paired_radial_defect=e,
                            signed_error_by_warm_parent=q.tolist()))
    return dict(x=x,target_q_lower_edges=list(PARENTS),parent_edges=edw.tolist(),moments=summary,
                maximum_change_in_untreated_cells=float(drift),
                scope='Positive targeted frozen quadrature at fixed occupations. Lower estimated source error does not establish pressure convergence or a better evolved model.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('receipt',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.state),read(a.receipt)),indent=2,allow_nan=False))
