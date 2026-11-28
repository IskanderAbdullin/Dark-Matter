"""Attribute a cold aperture change along a declared frozen-state replacement path."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Species
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import spacelike_parent_minimum


K=.0015


def coefficients(ch,m,omega,res,q,callbacks=None):
    Es=np.hypot(m,q)
    out=[]
    for a,b,c in ((4,6,12),(8,12,20)):
        out.append(angle(ch,K,omega,res,Es,'T',a,(b,c),q,occupation_callbacks=callbacks)[0])
    defect=float(max(abs(out[0]/out[1]-1)))
    if min(out[-1])<=0 or defect>2e-6:
        raise ArithmeticError(('Static counterfactual underresolved',q,m,omega,defect))
    return out[-1],defect


def run(cards):
    results=[]
    for g in sorted({'fine','third'}.intersection(c['grid'] for c in cards)):
        archive=load(g)[0]
        x=archive['x_end']
        base=[Species(mm,archive[f'endpoint_{name}_momentum_edges'],
                      archive[f'endpoint_{name}_occupations'])
              for mm,name in ((x,'parent'),(x*archive['mass_ratio'],'fermion'))]
        native_m=archive['checkpoints'][-1]['scalar_pole_over_T']
        _,w0,R0=root(lambda z:sum(np.array(c.field(z)) for c in base),K)
        for card in cards:
            if card['grid']!=g:
                continue
            if card['x_end']!=x:
                raise ArithmeticError('Compare only matched final epochs')
            d=CellChart(archive,SPECS[g],card['subcells_per_archived_cell'],reference_x=card['reference_x'])
            y=np.asarray(card['endpoint_occupation'])
            if card['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,):
                raise ArithmeticError('Invalid endpoint grid scope')
            m=d.state(x,y)[0]
            fit=[]
            ch=[]
            callbacks=[]
            for mass,e,n in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
                p,F,_,_,a,b=continuum_moments(e,mass,n,36)
                fit.append((p,F,mass))
                ch.append(piecewise_logit(mass,e,a,b,n,1))
                callbacks.append(model_callback(mass,e,a,b))
            _,w,R=root(lambda z:sum(charge_field(*k,z) for k in fit),K)
            for q in (3.5e-6,8.2e-6):
                expr=[(base,native_m,w0,R0,None),(base,m,w0,R0,None),
                      (base,m,w,R,None),(ch,m,w,R,callbacks)]
                data=[];err=[];geom=[]
                for b,mu,op,re,cl in expr:
                    A,F=coefficients(b,mu,op,re,q,cl)
                    data.append(A)
                    err.append(F)
                    t=np.hypot(mu,q)-op
                    V=np.hypot(K,q)
                    geom.append(spacelike_parent_minimum(x,x*d.ratio,t,V,(t-V)*(t+V)))
                z=np.asarray(data)
                ratio=z[1:]/z[:-1]
                if np.max(abs(np.prod(ratio,axis=0)/(z[-1]/z[0])-1))>3e-15:
                    raise ArithmeticError('Sequential path factors did not close')
                results.append(dict(grid=g,q=q,subcells=card['subcells_per_archived_cell'],
                                    original_and_new_scalar_poles=[native_m,m],
                                    original_and_new_photon_energy=[w0,w],
                                    raw_affine_coefficients=z.tolist(),mesh_discrepancies=err,
                                    radial_parent_threshold_for_cosine_zero=geom,
                                    replacement_order=['new_scalar_mass','new_photon_pole','new_continuum_charges'],
                                    per_stage_coefficient_factors=ratio.tolist(),
                                    total_factor=(z[-1]/z[0]).tolist()))
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cards',nargs='+',type=Path)
    o=p.parse_args()
    print(json.dumps(run([json.loads(x.read_text()) for x in o.cards]),indent=2,allow_nan=False))
