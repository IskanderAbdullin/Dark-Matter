"""Compare the full or mirror-reduced azimuth on physical retained polar branches."""
import json
from pathlib import Path
from time import process_time

import numpy as np

from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole,radial
from finite_mass_qed_timelike_pole_angular_reference import angle


def run():
    out=[]
    for g in ('fine','third'):
        d,ch,_=load(g)
        m=d['checkpoints'][-1]['scalar_pole_over_T']
        for name,k,q,caz,full in (('T',.0015,3.5e-6,0.,True),('T',.0015,.0009,-.61,False),
                                 ('T',.0015,.0021,.23,False),('T',.0015,.007,-.22,True),
                                 ('L',.0012,3.5e-6,.53,True),('L',.00256,.0012,-.37,False)):
            w,R=pole(ch,name,k)
            E=np.hypot(m,q)
            r=[]
            for mirror in (False,True):
                t=process_time()
                o,ok,err=radial(ch,k,w,R,E,caz,name,8,24,full,qscalar=q,
                       allow_future_spacelike=name=='T',mirror_azimuth=mirror)
                r.append((o,ok,err,process_time()-t))
            error=np.max(abs(r[0][0]-r[1][0])/max(np.max(abs(r[0][0])),1e-50))
            if error>2e-12 or r[0][1]!=r[1][1] or abs(r[0][2]-r[1][2])>2e-11:
                raise ArithmeticError(('Mirror symmetry or vertex replay disagrees',g,name,k,q,error))
            out.append(dict(grid=g,mode=name,k=k,q=q,cosine=caz,forward_inverse_signed=r[0][0].tolist(),
                            nonempty=r[0][1],relative_defect=float(error),symmetric_cpu_speedup=r[0][3]/r[1][3]))
    cases=[]
    for g,k,q in [('fine',.0015,3.5e-6),('third',.0015,.007)]:
        a,ch,_=load(g)
        m=a['checkpoints'][-1]['scalar_pole_over_T']
        w,R=pole(ch,'T',k)
        e=np.hypot(m,q)
        values=[]
        for mirror in (False,True):
            u,c,n,err=angle(ch,k,w,R,e,'T',4,(5,18),q,allow_future_spacelike=True,mirror_azimuth=mirror)
            values.append((u,c,n,err))
        a,b=values
        defect=float(np.max(abs(a[0]-b[0]))/max(np.max(abs(a[0])),1e-50))
        if defect>2e-12 or not np.array_equal(a[1],b[1]) or a[2:]!=b[2:]:
            raise ArithmeticError(('Full polar integration changed the mirror check',g,q,defect))
        cases.append(dict(grid=g,q=q,angular_panels=a[2],relative_defect=defect))
    if sum(x['nonempty'] for x in out)<7 or not any(x['mode']=='L' and x['nonempty'] for x in out):
        raise ArithmeticError('Too few active pole and support tests')
    result=continuum()
    return dict(radial=out,complete_cosine=cases,evolved_continuum=result,
                scope='Same discrete azimuth quadrature by exact mirror pairing; not a new angular refinement.')


def continuum():
    from early_pole_cell_endpoint_transverse_audit import (
        charge_field,continuum_moments,model_callback,piecewise_logit,root)
    from early_pole_endpoint_scalar_refined_pole_points import read
    from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
    from early_pole_endpoint_stencil_remap_pilot import SPECS,load as archive

    path=Path(__file__).parent/'results'/'early-pole-fine-charge2-scalar2-x05-state-2026-11-28.json.gz'
    src=read(path)
    x=src['x_end']
    d=ScalarSubshellChart(archive('fine')[0],SPECS['fine'],2,2,src['cut'])
    y=np.asarray(src['endpoint_occupation'])
    ch=[]; cb=[]; fit=[]
    for mass,edges,N in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,F,_,_,a,b=continuum_moments(edges,mass,N,36)
        fit.append((p,F,mass))
        ch.append(piecewise_logit(mass,edges,a,b,N,1))
        cb.append(model_callback(mass,edges,a,b))
    k=.0015
    m=d.state(x,y)[0]
    _,w,R=root(lambda z:sum(charge_field(*c,z,36) for c in fit),k)
    q=4.0e-6
    u,v=[angle(ch,k,w,R,np.hypot(m,q),'T',4,(6,18),q,
               occupation_callbacks=cb,mirror_azimuth=j)[0] for j in (False,True)]
    error=float(np.max(abs(u/v-1)))
    if error>2e-12:
        raise ArithmeticError(('Evolved continuum angular mirror did not match',error))
    return dict(q=q,relative_defect=error,affine=u.tolist())


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
