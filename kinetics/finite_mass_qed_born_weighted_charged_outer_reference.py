"""Signed Born weighting of the bare positive-frequency charged corrections."""
import gzip,json
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
from scipy.special import expit,roots_legendre
from early_pole_decay_mode_clock_reference import Q, intervals
from finite_mass_qed_charged_closed_radial_pv_reference import full_radial

D=np.longdouble


def compare(grid,d):
    M=d['x_end'];m=d['mass_ratio']*M
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,s,m,Q)
    p=np.hypot(d['endpoint_parent_momentum_edges'],M)
    b=np.hypot(d['endpoint_fermion_momentum_edges'],m)+Es
    fp=np.asarray(d['endpoint_parent_occupations'],dtype=D)
    fb=np.asarray(d['endpoint_fermion_occupations'],dtype=D)
    fs=D(d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1])
    frac=1/(1+fs)
    allnodes=sorted({lo,hi,*[x for x in (*p,*b) if lo<x<hi]})
    def weights(E):
        P=fp[np.searchsorted(p,E)-1]
        B=fb[np.searchsorted(b,E)-1]
        e=P*(1-B);a=B*(1-P)*(1-frac);net=(P-B)+frac*B*(1-P)
        if abs(e-a-net)>D('2e-18'):
            raise ArithmeticError('Signed weight failed its independently reconstructed detailed balance algebra.')
        return np.array([e,a,net])
    thermal=0.
    for E in np.linspace(lo,hi,11):
        ap,ab=expit(-E),expit(-(E-Es))
        a=ap*(1-ab)
        bth=ab*(1-ap)*np.exp(-Es)
        thermal=max(thermal,abs(a/bth-1))
    if thermal>2e-14:
        raise ArithmeticError('The equilibrium multiplicative virtual correction did not preserve equal forward and inverse weights.')

    @lru_cache(None)
    def bare(E,cut):
        R=full_radial(d,E,cut,lo,hi,False,local_aperture=True)
        return [R[k]['value'] for k in ('step','linear_energy')]

    rows=[]
    for tol in (2e-3,5e-4):
        a,c=lo+tol*(hi-lo),hi-tol*(hi-lo)
        splits=sorted({a,c,*[r for r in allnodes if a<r<c]})
        wave=[]
        for order in (5,9):
            nodes,w=roots_legendre(order)
            y=np.pi*(nodes+1)/2
            left=[]
            total=np.zeros((3,4),dtype=D)
            shape=np.zeros(3,dtype=D)
            for L,R in zip(splits[:-1],splits[1:],strict=True):
                E=L+(R-L)*np.sin(y/2)**2
                W=w*(R-L)*np.pi/4*np.sin(y)
                weight=weights((L+R)/2)
                A=np.asarray([bare(e,'parent')+bare(e,'daughter') for e in E])
                total+=np.outer(weight,W@A)
                shape+=weight*D(R-L)
                left.append(dict(energy=[L,R],archived_weights=[float(x) for x in weight]))
            if max(abs(total[0]-total[1]-total[2]))>2e-13 or min(shape)<=0:
                raise ArithmeticError('The positive, inverse and signed weighted projections did not close.')
            wave.append(dict(order=order,raw_weights=shape.astype(float).tolist(),weighted_averages=(total/shape[:,None]).astype(float).tolist(),cells=left))
        defect=float(np.max(np.abs(np.asarray(wave[1]['weighted_averages'])-wave[0]['weighted_averages'])))
        if defect>.04:
            raise ArithmeticError(f'Outer Born profile failed to resolve: {defect}.')
        rows.append(dict(outer_guard_fraction=tol,maximum_outer_resolution_defect=defect,controls=wave))
        print(f'finished {grid} external guard {tol:g}',file=sys.stderr,flush=True)
    return dict(grid=grid,scalar_occupation=float(fs),thermal_balance_relative_defect=thermal,rows=rows)


def run():
    results=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        results.append(compare(grid,d))
    return dict(scope='Bare positive-frequency charged correction on vacuum Born support, external weights normalized by 1+f_scalar. No physical cut coefficients or full rate.',
                weight_order=['forward','inverse','net'],
                vertex_columns=['parent_step','parent_linear_energy','daughter_step','daughter_linear_energy'],
                grids=results)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
