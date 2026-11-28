"""Test entropy-maximizing closure within a three-subcell charged bin."""
import json
from math import fsum

import numpy as np
from scipy.special import xlogy

from early_pole_adaptive_face_extents_pilot import once
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_bounded_face_extents_pilot import energy_difference,hold_trace


def project(f,w,E):
    f,w,E=(np.asarray(v).reshape(-1,3) for v in (f,w,E))
    v=np.column_stack([E[:,1]-E[:,2],E[:,2]-E[:,0],E[:,0]-E[:,1]])
    v/=np.max(abs(v),axis=1)[:,None]
    a=v/w
    lo=np.max(np.minimum(-f/a,(1-f)/a),axis=1)
    hi=np.min(np.maximum(-f/a,(1-f)/a),axis=1)
    if np.any(lo>hi):
        raise ArithmeticError('Infeasible Pauli/energy closure')
    for _ in range(55):
        m=(lo+hi)/2.
        p=np.clip(f+m[:,None]*a,np.nextafter(0.,1.),np.nextafter(1.,0.))
        gradient=np.sum(v*(np.log1p(-p)-np.log(p)),axis=1)
        lo=np.where(gradient>0,m,lo)
        hi=np.where(gradient>0,hi,m)
    alpha=(lo+hi)/2.
    z=f+alpha[:,None]*a
    for _ in range(5):
        bad=np.any((z<0.)|(z>1.),axis=1)
        if not np.any(bad):
            break
        alpha[bad]=np.nextafter(alpha[bad],0.)
        z=f+alpha[:,None]*a
    if z.min()<0. or z.max()>1.:
        raise ArithmeticError('Moment projection left its feasible interval')
    D=w*(z-f)
    if np.max(abs(D.sum(1)))>7e-15 or np.max(abs((D*E).sum(1)))>2e-14:
        raise ArithmeticError('Moment projection lost a cell invariant')
    entropy=lambda x:-xlogy(x,x)-xlogy(1-x,1-x)
    gain=fsum((w*(entropy(z)-entropy(f))).ravel())
    if gain<-2e-15:
        raise ArithmeticError('The constrained entropy closure decreased entropy')
    return z.ravel(),gain,fsum(abs(D).ravel())


def project_state(d,x,y):
    parts=np.split(y,(d.ns,d.iP))
    f,S,A=project(parts[1],d.d.wf,d.daughter_energies(x))
    p,T,B=project(parts[2],d.pw,d.parent_energies(x))
    return np.r_[parts[0],f,p],S+T,A+B


def run():
    ans=[]
    for grid,shape in SPECS.items():
        arc=load(grid)[0]
        x0=arc['x_end']
        for eta in (1/4096,1/32768):
            d=MassiveFermion(*shape,arc['mass_ratio'],arc['elastic_rate'],True)
            charge_mesh(d,eta,eta)
            y=np.concatenate([remap(a,arc['endpoint_'+k+'_occupations'],b)[0]
                              for a,b,k in zip(archives(arc),parts(d),LABELS,strict=True)])
            test=project_state(d,x0,y)[0]
            if np.max(abs(test-y))>5e-7:
                raise ArithmeticError('Initial constant archive not recovered by entropy closure')
            for n in (128,256):
                span=5e-5
                u=y.copy()
                es=[]
                mv=[]
                trace=[]
                h=span/n
                for j in range(n):
                    v,_,tr,_=once(d,x0,u,j*h,h)
                    # Constrained redistribution at the end of each
                    # completed bounded time step. Not a physical event.
                    xx=x0*np.exp((j+1)*h)
                    u,s,a=project_state(d,xx,v)
                    es.append(s)
                    mv.append(a)
                    trace.append(tr)
                diff=d.w*(u-y)
                nrg=energy_difference(d,x0,y,x0*np.exp(span),u)-fsum(trace)
                dn=[fsum(s) for s in np.split(diff,(d.ns,d.iP))]
                if max(abs(dn[0]-dn[1]),abs(dn[0]+dn[2]),abs(nrg))>1e-12:
                    raise ArithmeticError('The redistribution time pilot did not close')
                ans.append(dict(grid=grid,eta=eta,steps=n,paired_numbers=dn,
                                relative_face_motion_by_scalar_increment=fsum(mv)/max(abs(dn[0]),1e-50),
                                artificial_entropy_increase=fsum(es),dark_energy_residual=nrg))
    return ans


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
