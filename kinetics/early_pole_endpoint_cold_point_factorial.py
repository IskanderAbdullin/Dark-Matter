"""Separate gap, photon pole and charged closure in the cold-band's interior."""
import argparse
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from finite_mass_qed_timelike_pole_angular_reference import angle

K=.0015
FACTORS=('scalar_gap','photon_pole_including_residue','charged_occupation_closure')
POINTS=(2.5e-6,3.5e-6,4.7e-6,8.2e-6)


def prepare(c,d):
    if (c['grid'],c['x_end'],c['reference_x'],c['subcells_per_archived_cell'])!=('fine',.05,X,d):
        raise ArithmeticError('Need ordered fine d1 and d2 final states')
    a=load('fine')[0]
    x=c['x_end']
    chart=CellChart(a,SPECS['fine'],d,reference_x=X)
    y=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[chart.ns,chart.nf,chart.np] or y.shape!=(chart.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[chart.ns:])>1:
        raise ArithmeticError('Invalid state')
    pars=[]
    callbacks=[]
    ch=[]
    for mass,edges,n in ((x,chart.p_edges,chart.wp*y[chart.iP:]),(x*chart.ratio,chart.d_edges,chart.wd*y[chart.ns:chart.iP])):
        p,f,_,_,aa,b=continuum_moments(edges,mass,n,36)
        pars.append((p,f,mass))
        callbacks.append(model_callback(mass,edges,aa,b))
        ch.append(piecewise_logit(mass,edges,aa,b,n,1))
    _,omega,R=root(lambda u:sum(charge_field(*p,u,36) for p in pars),K)
    return chart.state(x,y)[0],(omega,R),(ch,callbacks),chart.ds.ed,y[:chart.ns]


def factors(values):
    log={k:np.log(v) for k,v in values.items()}
    if set(log)!=set(product((0,1),repeat=3)) or not np.all(np.isfinite(list(log.values()))):
        raise ArithmeticError('Incomplete positive factorial box')
    path=[]
    for order in permutations(range(3)):
        bit=[0,0,0]
        ans=np.empty((3,len(log[(0,0,0)])))
        for i in order:
            before=tuple(bit)
            bit[i]=1
            ans[i]=log[tuple(bit)]-log[before]
        path.append(ans)
    mat=np.stack(path)
    expectation=log[(1,1,1)]-log[(0,0,0)]
    if np.max(abs(mat.sum(axis=1)-expectation))>5e-14:
        raise ArithmeticError('The individual replacement paths fail to close')
    return np.mean(mat,axis=0),np.stack((np.min(mat,axis=0),np.max(mat,axis=0)))


def one(opts,q):
    def coeff(key,mesh):
        im,ip,ic=key
        m=opts[im][0]
        w,R=opts[ip][1]
        ch,cb=opts[ic][2]
        es=np.hypot(m,q)
        a,b,c=mesh
        u=angle(ch,K,w,R,es,'T',a,(b,c),q,occupation_callbacks=cb,mirror_azimuth=True)[0]
        if np.min(u)<=0 or abs(u[0]-u[1]-u[2])>4e-13*max(u):
            raise ArithmeticError('Invalid gross affine coefficients at a factorial corner')
        return u
    result={}
    rows=[]
    for key in product((0,1),repeat=3):
        a,b=[coeff(key,m) for m in ((6,10,18),(10,16,24))]
        e=np.max(abs(a/b-1))
        if e>2e-6:
            raise ArithmeticError(('Unresolved local quadrature',q,key,e))
        E0=np.hypot(opts[key[0]][0],q)
        c=b[:2]*q*q/E0
        result[key]=c
        rows.append(dict(bits=key,bare_forward_inverse=c.tolist(),mesh_defect=float(e)))
    for o in opts:
        if not o[3][0]<q<o[3][-1]:
            raise ArithmeticError('q outside the recorded scalar shell bounds')
    u,v=[o[4][np.searchsorted(o[3],q,side='right')-1] for o in opts]
    if min(u,v)<=0:
        raise ArithmeticError('No positive scalar occupations to compare')
    log_effect,path=factors(result)
    scalar=np.log(np.array([(1+v)/(1+u),v/u]))
    total=np.log(result[(1,1,1)]/result[(0,0,0)])+scalar
    if np.max(abs(log_effect.sum(axis=0)+scalar-total))>5e-14:
        raise ArithmeticError('The Bose/inverse scalar factor does not factorize')
    return dict(q=q,rows=rows,individual_path_log_factor_min_max=path.tolist(),
                mean_log_factors_by_bare_chart=log_effect.tolist(),scalar_occupations=[u,v],
                scalar_population_log_factor=scalar.tolist(),fully_weighted_d2_over_d1=np.exp(total).tolist())


def run(a,b):
    opts=[prepare(r,d) for r,d in ((a,1),(b,2))]
    if not np.array_equal(opts[0][3],opts[1][3]):
        raise ArithmeticError('Scalar shell edges changed in a charge-only comparison')
    return dict(k=K,x=.05,chart_factors=FACTORS,columns=('gross_forward','inverse'),
                chart_masses=[z[0] for z in opts],chart_photon_energy_residue=[z[1] for z in opts],
                points=[one(opts,q) for q in POINTS],
                scope='Four held-q replacement diagnostics. Arbitrary crossed states are not solutions, and point results are not q-integrated attribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('d1',type=Path)
    p.add_argument('d2',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.d1),read(a.d2)),indent=2,allow_nan=False))
