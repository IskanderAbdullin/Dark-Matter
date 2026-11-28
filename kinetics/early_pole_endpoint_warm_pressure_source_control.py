"""Separate frozen warm grid and evolved-state contrasts in an instantaneous pressure source."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def cell_weights(q,m,cut=None,order=32):
    a,b=q[:-1],q[1:]
    lo=a if cut is None else np.minimum(np.maximum(a,cut),b)
    h=b-lo
    t,w=roots_legendre(order)
    q0=lo[:,None]+h[:,None]*(t+1)/2
    total=(b-a)*(b*b+a*b+a*a)/3
    integral=(h/2)*((q0**4/np.hypot(q0,m))@w)
    return integral/total


def local(d,y,x,m,H,n):
    mesh,project=d.chart(x,m,n)
    vals,A,B=event(mesh,y,d)
    idx=mesh[3]
    p=d.energies(x,m)
    energy=fsum(vals*np.concatenate(p))
    number=[fsum(k) for k in np.split(vals,(d.ns,d.iP))]
    if max(abs(energy),abs(number[0]-number[1]),abs(number[0]+number[2]))>3e-13 or project>3e-14:
        raise ArithmeticError('Lost conservation in a frozen collision source')
    P,W=[cell_weights(d.ds.ed,m,c,64) for c in (None,.01)]
    for ref,cut in ((P,None),(W,.01)):
        coarse=cell_weights(d.ds.ed,m,cut,32)
        if np.max(abs(coarse-ref)/np.maximum(abs(ref),1e-105))>2e-13:
            raise ArithmeticError('The within-cell weights did not resolve')
    weights=(np.ones(d.ns),P,W)
    Y=np.array([[H*fsum(z[idx]*q) for q in (A,B,A-B)] for z in weights])
    if np.max(abs(Y[:,0]-Y[:,1]-Y[:,2])/np.maximum(np.max(abs(Y),axis=1),1e-99))>4e-14:
        raise ArithmeticError('Gross and signed source receipts do not close')
    return Y,project


def run(base,warm,x):
    (a,u),(b,w)=unpack(base,False,x),unpack(warm,True,x)
    if any(not np.array_equal(getattr(a,k),getattr(b,k)) for k in ('qd','qp','wd','wp')):
        raise ArithmeticError('The frozen comparison must keep the charge chart fixed')
    idx=np.searchsorted(a.ds.ed,b.ds.p)-1
    v=np.r_[u[idx],u[a.ns:]]
    if abs(fsum(a.ds.w*u[:a.ns])/fsum(b.ds.w*v[:b.ns])-1)>3e-14:
        raise ArithmeticError('The step-function lift lost scalar number')
    m,_,_,H=a.state(x,u)
    out=[]
    for d,y in ((a,u),(b,v),(b,w)):
        (lo,L),(hi,U)=[local(d,y,x,m,H,n) for n in (4,8)]
        error=np.max(abs(lo-hi)/np.maximum(abs(hi),1e-99))
        if error>1e-6:
            raise ArithmeticError(('The paired local quadratures have not resolved',error))
        out.append(dict(gross_forward_inverse_net=hi.tolist(),relative_mesh_defect=float(error),
                        eventwise_energy_projection=U))
    I,J,K=[np.array(o['gross_forward_inverse_net']) for o in out]
    return dict(x=x,held_mass=float(m),held_expansion_normalization=float(H),
                observables=['total_scalar_number','reconstructed_pressure_numerator','pressure_at_q_greater_equal_point01'],
                cases=['S4_original','W2_with_lifted_S4_step','W2_actual_state'],
                records=out,frozen_chart_relative_net_change=(J[:,2]/I[:,2]-1).tolist(),
                evolved_state_given_refined_chart_relative_net_change=(K[:,2]/J[:,2]-1).tolist(),
                combined_relative_net_change=(K[:,2]/I[:,2]-1).tolist(),
                scope='Instantaneous collision source of reconstructed moments at held mass and clock, not their total time derivatives.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('S4',type=Path)
    p.add_argument('W2',type=Path)
    p.add_argument('--x',type=float,choices=(.01,.02,.04,.05),default=.05)
    a=p.parse_args()
    print(json.dumps(run(read(a.S4),read(a.W2),a.x),indent=2,allow_nan=False))
