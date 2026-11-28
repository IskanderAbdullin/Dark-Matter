"""Measure retained parent-energy localization at fixed pole and cold scalar momenta."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import radial


EDGES=(0.,16.,18.,20.,22.,24.,26.,28.,30.,40.)


def measure(c):
    g=c['grid']
    x=c['x_end']
    if x!=.05 or c['subcells_per_archived_cell']!=4:
        raise ArithmeticError('Use a four-charge-subcell final state')
    d=CellChart(load(g)[0],SPECS[g],4,reference_x=c['reference_x'])
    y=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid final state')
    m=d.state(x,y)[0]
    fit=[]; ch=[]; cb=[]
    for mass,e,n in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,F,_,_,a,b=continuum_moments(e,mass,n,36)
        fit.append((p,F,mass))
        ch.append(piecewise_logit(mass,e,a,b,n,1))
        cb.append(model_callback(mass,e,a,b))
    k=.0015
    _,w,R=root(lambda beta:sum(charge_field(*a,beta,36) for a in fit),k)
    out=[]
    kw=dict(occupation_callbacks=cb,mirror_azimuth=True)
    for q in (3.5e-6,1.5e-5):
        E=np.hypot(m,q)
        for cosine in (None,.9):
            def sample(mesh,edges):
                n,e,az=mesh
                if cosine is None:
                    return angle(ch,k,w,R,E,'T',n,(e,az),q,parent_energy_edges=edges,**kw)[0]
                return radial(ch,k,w,R,E,cosine,'T',e,az,qscalar=q,parent_energy_edges=edges,**kw)[0]
            a,b=[sample(order,EDGES) for order in ((4,6,12),(8,10,20))]
            total=b.sum(axis=0)
            ref=sample((8,10,20),None)
            # Allow exact zero bands and assess tiny pieces against their column total.
            f=lambda p,q: float(np.max(abs(p-q)/np.maximum(abs(total),1e-80)))
            drift=f(a,b)
            closure=f(total,ref)
            if max(drift,closure)>2e-6 or total[0]<=0 or total[1]<0:
                raise ArithmeticError(('Unresolved energy tail partition',g,q,cosine,drift,closure))
            share=[(b[:,j]/total[j]).tolist() if total[j]>0 else None for j in (0,1)]
            above=[sum(u[EDGES.index(28.):]) if u is not None else None for u in share]
            if (share[1] is None and (np.any(a[:,1]!=0) or np.any(b[:,1]!=0) or ref[1]!=0)):
                raise ArithmeticError('An unresolved inverse is not certified as absent')
            out.append(dict(q=q,cosine=cosine,group_affine=b.tolist(),
                            fractions_of_forward_inverse=share,above_28_fractions=above,
                            mesh_absolute_fraction_drift=drift,unsplit_affine_relative_closure=closure))
    return dict(grid=g,k=k,state_scalar_mass=m,omega=w,parent_energy_edges=EDGES,points=out,
                scope='Retained parent-energy localization only, no extrapolation beyond archived occupations.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',type=Path,nargs='+')
    a=p.parse_args()
    print(json.dumps([measure(read(c)) for c in a.states],indent=2,allow_nan=False))
