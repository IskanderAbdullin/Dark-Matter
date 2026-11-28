"""Pick per-cell warm orders at x=.01 and test unchanged rules at x=.05."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_positive_gauss_ladder import volume_gauss
from early_pole_endpoint_warm_source_cells import local


def choose(loss,cap,orders):
    loss=np.asarray(loss,dtype=float)
    if loss.ndim!=2 or not loss.size or loss.shape[1]!=len(orders) or not np.all(np.isfinite(loss)) or np.min(loss)<0:
        raise ArithmeticError('Need nonnegative finite cell losses')
    states={0:(0.,())}
    for row in loss:
        t={}
        for n,(m,v) in states.items():
            for p,c in zip(orders,row,strict=True):
                j=n+p
                if j>cap:
                    continue
                promise=(m+float(c),(*v,p))
                if j not in t or promise<t[j]:
                    t[j]=promise
        states=t
    if not states:
        raise ArithmeticError('Too few quadrature nodes under the selected budget')
    n,z=min(states.items(),key=lambda kv:(kv[1][0],kv[0],kv[1][1]))
    return n,z


def table(meta,comparison):
    if meta['x']!=comparison['x'] or meta['warm_native_edges']!=comparison['archived_warm_edges']:
        raise ArithmeticError('Frozen finite source tables do not share a state')
    cases={r['rule']:np.asarray(r['warm_cell_forward_inverse_net']) for r in meta['cases']}
    typical=comparison['comparison_columns']
    N=meta['ref_net_pressure']
    if typical!=['forward','inverse','net'] or N!=comparison['reference_pressure_rate'] or N<=0.:
        raise ArithmeticError('Frozen rate normalization changed')
    levels={2:cases['Gauss2']/N}
    for row in comparison['changes']:
        a,b=row['orders']
        if b!=a+1 or b in levels or a not in levels:
            raise ArithmeticError('Need every adjacent direct order comparison once')
        levels[b]=levels[a]+np.asarray(row['all_warm_parent_changes'])
    V=cases['Volume64']/N
    for k,row in zip(range(2,9),comparison['charts'],strict=True):
        if row['order']!=k or abs(fsum(abs(levels[k][:,2]-V[:,2]))-row['l1_signed_source_residual'])>3e-12:
            raise ArithmeticError('Recovered local differences do not close the published whole totals')
    return levels,V


class Mixed(ScalarSubshellChart):
    def __init__(self,orders):
        super().__init__(load('fine')[0],SPECS['fine'],2,4,1e-5)
        a,b=self.ds.ed[:-1],self.ds.ed[1:]
        self.old_count=self.ns
        warm=a>=.01
        self.repeats=np.ones(self.ns,dtype=int)
        self.repeats[warm]=orders
        Q=[]
        W=[]
        for i,p in enumerate(self.repeats):
            if p==1:
                Q.append([self.ds.p[i]])
                W.append([self.ds.w[i]])
            else:
                q,w,_=volume_gauss(a[i:i+1],b[i:i+1],p)
                Q.append(q[0])
                W.append(w[0]*(self.ds.w[i]/w.sum()))
        self.ds.p=np.concatenate(Q)
        self.ds.w=np.concatenate(W)
        self.ds.n=self.ns=len(self.ds.p)
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]

    def lift(self,y):
        return np.r_[np.repeat(y[:self.old_count],self.repeats),y[self.old_count:]]


def evaluate(card,ref,comp,choice):
    x=card['x_end']
    a,y=unpack(card,False,x)
    m,_,_,H=a.state(x,y)
    q=Mixed(choice)
    u=q.lift(y)
    warm=a.ds.ed[:-1]>=ref['warm_floor']
    if (x,m,H)!=(ref['x'],ref['mass_held'],ref['clock_held']) or any(not np.array_equal(getattr(a,k),getattr(q,k)) for k in ('wp','wd','qp','qd')):
        raise ArithmeticError('Changed the independently stored optics, clock or charge charts')
    f=np.add.reduceat(q.ds.w*u[:q.ns],np.r_[0,np.cumsum(q.repeats)[:-1]])
    if np.max(abs(f-a.ds.w*y[:a.ns])/np.maximum(a.ds.w*y[:a.ns],1e-105))>3e-14:
        raise ArithmeticError('Mixed lifting changed native scalar density')
    P=np.where(q.ds.p>=ref['warm_floor'],q.ds.p**2/np.hypot(q.ds.p,m),0.)
    l,r=[local(q,u,P,x,m,H,a.ds.ed,n)[warm] for n in (4,8)]
    N=ref['ref_net_pressure']
    F,V=table(ref,comp)
    expected=np.array([F[o][k] for k,o in enumerate(choice)])
    v=r/N
    defect=float(np.max(abs(l-r)/np.maximum(r[:,0,None]+r[:,1,None],1e-100)))
    local_max=float(np.max(abs(v-expected)))
    if defect>3e-6 or local_max>3e-12 or q.ns!=a.ns-len(choice)+sum(choice):
        raise ArithmeticError('Fresh mixed source does not reproduce independent held per-cell assembly')
    err=v[:,2]-V[:,2]
    return dict(x=x,nodes=q.ns,reference_l1=fsum(abs(err)),signed_residual=fsum(err),
                maximum_radial_mesh_defect=defect,
                maximum_absolute_mixed_minus_independent_per_cell=local_max,
                individual_cell_net_residual=err.tolist())


def run(states,reference,comparisons):
    if [z['x_end'] for z in states]!=[.01,.05] or reference[0]['warm_native_edges']!=reference[1]['warm_native_edges']:
        raise ArithmeticError('Require one prespecified training and one separate held check endpoint')
    orders=(4,5,6,7,8)
    k,V=table(reference[0],comparisons[0])
    e=abs(np.array([k[p][:,2]-V[:,2] for p in orders])).T
    a,_=unpack(states[0],False,.01)
    cold=a.ns-len(e)
    choices=[]
    for budget in (219,255):
        n,(loss,path)=choose(e,budget-cold,orders)
        z=[evaluate(a,b,c,path) for a,b,c in zip(states,reference,comparisons,strict=True)]
        if abs(z[0]['reference_l1']-loss)>3e-12:
            raise ArithmeticError('Profile at the training state did not reproduce the selected loss')
        choices.append(dict(maximum_total_scalar_nodes=budget,applied_scalar_nodes=n+cold,
                            warm_orders_by_cell=path,paired_direct_results=z))
    return dict(training_x=.01,held_check_x=.05,orders=orders,warm_edges=reference[0]['warm_native_edges'],
                choices=choices,scope='Optimizes a frozen numerical source at one endpoint, then directly evaluates another fixed endpoint on the same baseline path. Not evolved adaptation or statistical generalization.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('states','reference','comparisons'):
        p.add_argument('--'+key,type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(run(*[[read(f) for f in getattr(a,k)] for k in ('states','reference','comparisons')]),indent=2,allow_nan=False))
