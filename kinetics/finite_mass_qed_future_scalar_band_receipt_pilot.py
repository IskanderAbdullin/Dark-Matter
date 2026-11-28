"""Compare complete opt-in future and cold band matrices without filling absent q."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


def read(path):
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',encoding='utf-8') as f:
        return [json.loads(s) for s in f]


def combine(path,grid):
    data=read(path)
    if not data:
        raise ArithmeticError(('Empty future-band receipt',path))
    total=data[0]['total']
    if len(data)!=total or [x['panel'] for x in data]!=list(range(total)):
        raise ArithmeticError(('Partially completed q intervals',path,len(data),total))
    spans=[x['limits'] for x in data]
    if (spans[0][0],spans[-1][1])!=(.004,.015) or any(
         a[1]!=b[0] for a,b in zip(spans[:-1],spans[1:],strict=True)):
        raise ArithmeticError(('Future raw panels contain holes',path))
    column=('mesh','band','k','omega','total')
    if any(x['grid']!=grid or any(x[key]!=data[0][key] for key in column) for x in data) or data[0]['k']!=.0015:
        raise ArithmeticError('Wrong archive/grid/mode in receipt')
    a=np.array([z['number_scalar_energy_and_centered_coefficients'] for z in data])
    if np.min(a[:,:,:2])<0 or not np.all(np.isfinite(a)):
        raise ArithmeticError('Negative gross process weight')
    defect=abs(a[:,:,0]-a[:,:,1]-a[:,:,2])/np.maximum(np.max(abs(a),axis=2),1e-50)
    if np.max(defect)>5e-13:
        raise ArithmeticError(('Affine decomposition did not close',np.max(defect)))
    return np.array([[fsum(a[:,i,j]) for j in range(3)] for i in range(3)]),spans,data[0]['omega']


def run(root,cold,pack,meshes=(0,1)):
    c={x['grid']:x for x in json.loads(cold.read_text())}
    out=[]
    for g in ('fine','third'):
        a,b=[combine(root/f'{g}-first-future-{i}.jsonl{pack}',g) for i in meshes]
        if a[1]!=b[1] or a[2]!=b[2] or b[2]!=c[g]['pole_omega']:
            raise ArithmeticError(('Refinement or cold band changed support/pole',g))
        x=b[0]
        old=np.asarray(c[g]['integrated_number_scalar_energy_and_excess_coefficients'])
        # Pair the full and centered scalar energy to infer the same
        # scalar rest mass on each independently assembled source column.
        rest=(x[1]-x[2])/x[0]
        oldrest=(old[1]-old[2])/old[0]
        if np.max(abs(rest-oldrest))/np.max(oldrest)>4e-13:
            raise ArithmeticError(('Inconsistent energy centering',rest,oldrest))
        e=abs(a[0]-x)/np.maximum(np.maximum(abs(a[0]),abs(x)),1e-45)
        if np.max(e)>3e-5:
            raise ArithmeticError(('Unresolved future-q mesh',g,e))
        if np.any(x[:2,1:]<=0):
            raise ArithmeticError(('No positive conditional threshold',g,x))
        theta=x[:2,1]/x[:2,2]
        ctheta=c[g]['number_and_scalar_energy_threshold']
        # Keep either band distinct: this is not a prediction of actual n.
        w=b[2]
        if max(w,c[g]['pole_omega'])>=.004:
            raise ArithmeticError('Future region no longer has positive t throughout')
        averaged_energies=x[1,:2]/x[0,:2]
        if min(averaged_energies)<=.004 or max(averaged_energies)>=.016:
            raise ArithmeticError('Gross scalar energies left the declared band')
        trials=1e-8,1e-5,*theta,float(np.mean(theta))
        fits=[]
        for n in trials:
            cold_flux=n*old[:2,2]-old[:2,1]
            delta_n,delta_s=n*x[:2,2]-x[:2,1]
            delta_gamma=-w*delta_n
            delta_charges=-delta_s-delta_gamma
            # Reconstruct charged transfer directly from the two extensive
            # source columns, not from the already combined scalar changes.
            raw_charges=n*(w*x[0,2]-x[1,2])-(w*x[0,1]-x[1,1])
            if abs(delta_charges-raw_charges)>2e-13*max(abs(delta_s),abs(delta_gamma),1e-50):
                raise ArithmeticError(('Charged gain or loss does not close',g,n,delta_charges,raw_charges))
            fits.append(dict(n=n,cold=cold_flux.tolist(),future=[delta_n,delta_s],
                             future_photon_energy=delta_gamma,future_charged_energy=delta_charges))
        out.append(dict(grid=g,future_panels=len(b[1]),coefficient_drift=e.tolist(),
                        coefficients=x.tolist(),number_and_scalar_energy_neutral=theta.tolist(),
                        ratio_of_future_to_cold_forward=(x[:2,0]/old[:2,0]).tolist(),
                        cold_number_energy_neutral=ctheta,mean_scalar_energies_forward_inverse=averaged_energies.tolist(),
                        comparative_sources=fits))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).parent/'results'
    p.add_argument('--root',type=Path,default=here)
    p.add_argument('--cold',type=Path,default=here/'finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json')
    p.add_argument('--raw',action='store_true')
    p.add_argument('--meshes',nargs=2,type=int,choices=(0,1,2),default=(0,1))
    o=p.parse_args()
    if not o.meshes[0]<o.meshes[1]:
        raise ValueError('Mesh refinement must increase quadrature order')
    print(json.dumps(run(o.root,o.cold,'' if o.raw else '.gz',o.meshes),indent=2,allow_nan=False))
