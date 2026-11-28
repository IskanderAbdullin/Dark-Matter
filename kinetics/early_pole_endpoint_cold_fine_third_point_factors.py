"""Separate observed fine/third D2 differences at declared cold q points."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_charge_cold_clock import optical
from early_pole_endpoint_charge_cold_clock_points import point,setup
from early_pole_endpoint_cold_non_nested_cdf_contrast import audit
from early_pole_endpoint_cold_point_factorial import FACTORS,POINTS,one
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(fine,third,states,archive,cdf):
    if cdf['grids']!=['fine','third'] or cdf['steps']!=2048 or cdf['charge_subcells']!=2:
        raise ArithmeticError('Expected independently validated non-nested whole cold bands')
    drafts=[]
    for g,pair,path,sha in zip(('fine','third'),(fine,third),states,cdf['states_sha256_uncompressed'],strict=True):
        _,sig,receipt=audit(pair,path,archive,g)
        if sig!=sha or receipt!=cdf['gates'][len(drafts)]:
            raise ArithmeticError('The point factors do not bind the already accepted whole bands')
        c,_=read(path)
        d=CellChart(load(g)[0],SPECS[g],2,reference_x=X)
        y=np.asarray(c['endpoint_occupation'])
        if (c['reference_x']!=X or c['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,)
            or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1):
            raise ArithmeticError('Wrong or nonpositive native optical state')
        m,w,res=optical(d,y,.05)
        a=receipt['full_pole_and_snapshot']
        if np.max(abs(np.asarray([m,w,res])/[a[k] for k in ('scalar_mass','omega','residue')]-1.))>2e-13:
            raise ArithmeticError('Optical endpoint fit does not reproduce the archived full-band pole')
        rig=setup(d,y,(m,w,res))
        ch,cb,*_=rig
        drafts.append((d,y,rig,(m,(w,res),(ch,cb),d.ds.ed,y[:d.ns])))
    if drafts[0][0].ratio!=drafts[1][0].ratio:
        raise ArithmeticError('Unable to isolate charts with different reaction masses')
    out=[]
    for q in POINTS:
        row=one([z[-1] for z in drafts],q)
        ratio=row.pop('fully_weighted_d2_over_d1')
        s,t=[point(d,y,r,q) for d,y,r,_ in drafts]
        native=t[0][:2]/s[0][:2]
        proof=float(np.max(abs(np.asarray(ratio)/native-1.)))
        if max(proof,s[1],t[1])>2e-9:
            raise ArithmeticError('Factorized diagonal did not reproduce independent fully weighted points')
        col=np.r_[row['mean_log_factors_by_bare_chart'],[row['scalar_population_log_factor']]]
        share=col/np.log(native)
        out.append(dict(**row,fully_weighted_third_over_fine=native.tolist(),
                        bare_and_explicit_scalar_log_shares=share.tolist(),
                        maximum_independent_point_diagonal_relative_defect=proof))
    return dict(x=.05,k=.0015,grids=['fine','third'],charge_subcells=2,fixed_final_steps=2048,
                states_sha256_uncompressed=cdf['states_sha256_uncompressed'],fields=['forward_number','inverse_number'],
                factor_order=[*FACTORS,'explicit_scalar_population'],points=out,
                scope='Independent point-only frozen numerical substitutions. The full-band CDF/overlap comparison remains separate; no causal or continuum inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','third','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    for k in ('archive','cdf'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive,json.loads(a.cdf.read_text())),indent=2,allow_nan=False))
