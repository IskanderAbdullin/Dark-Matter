"""Average the doubly conservative radial reconstruction under signed Born weights."""
from functools import lru_cache
import gzip,json
from pathlib import Path
import sys

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals
from finite_mass_qed_conservative_virtual_reconstruction_reference import Smoother, calc

D=np.longdouble


def compare(d,g,witness):
    M=d['x_end'];m=d['mass_ratio']*M
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,s,m,Q)
    edges={k:np.asarray(d[f'endpoint_{k}_momentum_edges']) for k in ('parent','fermion')}
    occ={k:np.asarray(d[f'endpoint_{k}_occupations']) for k in ('parent','fermion')}
    Eparent=np.hypot(M,edges['parent'])
    Echild=np.hypot(m,edges['fermion'])+Es
    fs=D(d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1])
    fplus=1/(1+fs)
    width=.05
    sp=Smoother(edges['parent'],occ['parent'],M,width)
    sb=Smoother(edges['fermion'],occ['fermion'],m,width)

    @lru_cache(None)
    def grab(E):
        pa=calc(d,E,'parent',lo,hi,width,False,sp,orders=(48,72))
        fb=calc(d,E,'daughter',lo,hi,width,False,sb,orders=(48,72))
        return [pa['bare_virtual'],fb['bare_virtual']],max(pa['radial_mesh_defect'],fb['radial_mesh_defect'])

    out=[]
    for guard in (.002,.0005):
        L,H=lo+guard*(hi-lo),hi-guard*(hi-lo)
        split=sorted({L,H,*[x for x in (*Eparent,*Echild,*sp.boundaries.astype(float),*(sb.boundaries.astype(float)+Es)) if L<x<H]})
        shaped=[]
        for n in (6,10):
            t,w=roots_legendre(n)
            u=np.pi*(t+1)/2
            x=np.sin(u/2)**2
            V=w*np.pi*np.sin(u)/4
            y=np.zeros((3,2),dtype=D);den=np.zeros(3,dtype=D)
            radial=0.
            for l,h in zip(split[:-1],split[1:],strict=True):
                c=(l+h)/2
                fp=D(occ['parent'][np.searchsorted(Eparent,c)-1])
                fb=D(occ['fermion'][np.searchsorted(Echild,c)-1])
                weights=np.array([fp*(1-fb),fb*(1-fp)*(1-fplus),fp-fb+fb*(1-fp)*fplus])
                pieces=[grab(float(v)) for v in l+(h-l)*x]
                vals=np.asarray([o[0] for o in pieces])
                radial=max(radial,max(o[1] for o in pieces))
                y+=(h-l)*np.outer(weights,V@vals)
                den+=(h-l)*weights
            if max(abs(y[0]-y[1]-y[2]))>2e-13:
                raise ArithmeticError('External common-matrix forward/inverse/net replay failed.')
            shaped.append(dict(order=n,coefficients=(y/den[:,None]).astype(float).tolist(),
                               raw_normalized_weights=den.astype(float).tolist(),worst_radial_mesh=radial))
        diff=np.asarray(shaped[-1]['coefficients'])-shaped[0]['coefficients']
        old=next(k for k in witness['rows'] if k['outer_guard_fraction']==guard)['controls'][-1]
        if max(abs(np.array(old['raw_weights'])-np.array(shaped[-1]['raw_normalized_weights'])))>2e-15:
            raise ArithmeticError('Signed external weights did not match the archive reference.')
        defect=np.max(abs(diff))
        if defect>2e-4:
            raise ArithmeticError(f'Conservative outer average did not converge: {defect}.')
        z=np.asarray(shaped[-1]['coefficients']).sum(axis=1)
        base=np.asarray(old['weighted_averages'])
        step=base[:,0]+base[:,2]
        linear=base[:,1]+base[:,3]
        row=dict(endpoint_guard_fraction=guard,outer_mesh_defect=float(defect),samples=shaped,
                 step_cuts_signed_total=step.tolist(),simple_linear_cuts_signed_total=linear.tolist(),
                 doubled_conservative_cuts_signed_total=z.tolist(),
                 relative_change_from_step=(z/step-1).tolist())
        print(g,guard,row['relative_change_from_step'],file=sys.stderr,flush=True)
        out.append(row)
    return dict(grid=g,request_transition_fraction=width,weight_order=('forward','inverse','net'),rows=out)


def run():
    root=Path(__file__).parent/'results'
    with open(root/'finite-mass-qed-born-weighted-charged-outer-2026-11-27.json') as f: prior=json.load(f)
    out=[]
    for g in ('fine','third'):
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        out.append(compare(d,g,next(x for x in prior['grids'] if x['grid']==g)))
    return dict(scope='Signed external Born occupations and internal Pauli-bounded continuous reconstruction preserving each cell number and energy. Isolated bare positive charged cuts only.',grids=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
