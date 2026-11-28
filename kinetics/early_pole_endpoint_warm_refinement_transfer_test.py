"""Test if cheap warm refinement choices transfer to independent archived times."""
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


def source(c,bundle,reference):
    x=c['x_end']
    a,y=unpack(c,False,x)
    m,_,_,H=a.state(x,y)
    ed=a.ds.ed
    i=np.searchsorted(ed,reference['warm_floor'])
    e=ed[i:]
    if (x,m,H)!=(reference['x'],reference['mass_held'],reference['clock_held']) or e.tolist()!=reference['warm_native_edges']:
        raise ArithmeticError('Frozen references are not matched')
    if reference['mesh_reference']['global_L1_fraction']>4e-5:
        raise ArithmeticError('Warm volume reference did not close')
    ref={r['rule']:np.asarray(r['warm_cell_forward_inverse_net']) for r in reference['cases']}
    R=ref['Volume64']
    T=fsum(R[:,2])
    if T<=0 or abs(T/reference['ref_net_pressure']-1)>2e-14:
        raise ArithmeticError('Use the positive accepted reference scale')
    N=a.ds.w*y[:a.ns]
    def calc(chart):
        z=chart.transfer_step(y)
        pos=np.r_[0,np.cumsum(chart.repeats)[:-1]]
        n=np.add.reduceat(chart.ds.w*z[:chart.ns],pos)
        if (np.max(abs(n-N)/np.maximum(N,1e-105))>3e-14 or not np.array_equal(chart.ds.ed,ed)
                or any(not np.array_equal(getattr(a,k),getattr(chart,k)) for k in ('wp','wd','qp','qd'))):
            raise ArithmeticError('Changed archived counts or charge grid')
        q=chart.ds.p
        P=np.where(q>=e[0],q*q/np.hypot(q,m),0.)
        lo,hi=[local(chart,z,P,x,m,H,ed,k)[i:] for k in (4,8)]
        g=np.max(abs(lo-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100))
        if g>3e-6:
            raise ArithmeticError('Radial mesh not resolved')
        return hi
    out=[calc(chart) for chart in (WarmGauss(bundle,4),WarmGauss(bundle,5),Targeted(bundle,e[:-1]))]
    def validate(k):
        keep=np.zeros(len(e)-1,dtype=bool)
        keep[k]=True
        want=np.where(keep[:,None],out[2],out[0])
        actual=calc(Targeted(bundle,e[:-1][keep]))
        return float(np.max(abs(want-actual)/np.maximum(np.max(abs(want),axis=1)[:,None],1e-100)))
    return dict(x=x,edges=e,reference=R[:,2],reference32=ref['Volume32'][:,2],scale=T,base=out[0][:,2],uniform5=out[1][:,2],split=out[2][:,2],pop=N[i:],validate=validate)


def stats(q,mask,tag):
    k=np.zeros_like(q['base'],dtype=bool)
    k[np.asarray(mask,dtype=int)]=True
    pred=np.where(k,q['split'],q['base'])
    if tag=='uniform_G5':
        pred=q['uniform5']
    D=(pred-q['reference'])/q['scale']
    low=(pred-q.get('reference32',q['reference']))/q['scale']
    return dict(rule=tag,L1_against_Volume32_in_Volume64_scale=fsum(abs(low)),refined_parent_edges=q['edges'][:-1][k].tolist(),signed=fsum(D),
                sum_abs=fsum(abs(D)),worst_abs=float(np.max(abs(D))),
                chart_nodes=233 if tag=='uniform_G5' else 211+4*int(np.count_nonzero(k)))


def run(cards,refs):
    g=load('fine')[0]
    Q=[source(c,g,r) for c,r in zip(cards,refs,strict=True)]
    if len(Q)!=4 or [x['x'] for x in Q]!=[.01,.02,.04,.05] or any(not np.array_equal(q['edges'],Q[0]['edges']) for q in Q):
        raise ArithmeticError('Need the four predeclared S4 benchmarks')
    a=Q[0]
    def top(v):
        return np.argsort(-v,kind='stable')[:2]
    order={'published_pair':np.searchsorted(a['edges'],PARENTS),
           'train01_split_discrepancy':top(abs(a['split']-a['base'])),
           'train01_G4_source_magnitude':top(abs(a['base'])),
           'train01_population':top(a['pop'])}
    results=[]
    for q in Q:
        o=[stats(q,[],'uniform_G4'),stats(q,[],'uniform_G5')]
        o.extend(stats(q,M,k) for k,M in order.items())
        err=q['validate'](order['train01_split_discrepancy'])
        if err>3e-13:
            raise ArithmeticError('Whole chart did not reproduce independent replacement composition')
        results.append(dict(x=q['x'],validation=q['x']!=.01,moments=o,assembled_vs_actual_selected_maximum_defect=err))
    return dict(training_x=.01,comparisons=results,
                scope='Fixed early-time ranking, no validation labels in the ranking. Replace only the two selected native cell source integrals. This does not evolve either chart.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    p.add_argument('--references',nargs=4,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run([read(p) for p in a.states],[read(p) for p in a.references]),indent=2,allow_nan=False))
