"""Bound the Hartree root of any frozen positive shape with reference cell moments."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_endpoint_scalar_energy_moment_envelope import envelopes
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import fit,scaled_energy,shell_rule
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def case(card,lift):
    g=card['grid']
    d=ScalarGauss2(load(g)[0],SPECS[g],card['charge_div'],card['scalar_div'],card['cut'],
                   gauss_floor=card.get('gauss_floor'))
    x=card['x_end']
    if lift and (g,x,card['scalar_div'],card['counts'][0])!=('fine',.05,8,d.old_ns):
        raise ArithmeticError('Use the specified matched S8 frozen histogram')
    if card['reference_x']!=X or (not lift and card.get('scalar_nodal_rule')!='gauss-volume-v1'):
        raise ArithmeticError('Unsupported Gauss input chart')
    y=d.transfer_step(np.asarray(card['endpoint_occupation'])) if lift else np.asarray(card['nodal_occupation'])
    if y.shape!=(d.n,) or min(y)<0 or not np.all(np.isfinite(y)) or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid positive nodal population')
    ref=d.state(x,y)[0]
    vol=d.ds.w*y[:d.ns]
    a,b=d.ds.ed[:-1],d.ds.ed[1:]
    V=np.r_[0,np.cumsum(d.scalar_repeats)]
    N=[]; f=[]; E=[]; Q=[]; W=[]; means=[]
    for i,(l,r) in enumerate(zip(a,b,strict=True)):
        ix=slice(V[i],V[i+1])
        n=fsum(vol[ix])
        if not n:
            continue
        s=fsum(vol[ix]*scaled_energy(d.ds.p[ix],l,r,ref))/n
        nr=fit(l,r,ref,n,s)
        q,eta,z=shell_rule(l,r,ref,128)
        P=np.exp(z+nr['log_occupation_at_lower_edge']+nr['tilt']*eta)
        if abs(fsum(P)/n-1)>5e-12:
            raise ArithmeticError('Reconstructed fixed shape changed the count')
        lo=l*l/(np.hypot(l,ref)+ref)
        hi=r*r/(np.hypot(r,ref)+ref)
        N.append(n)
        f.append(lo)
        E.append(hi)
        Q.append(q)
        W.append(P)
        # Keep the mean in excess-energy coordinates, stable for cold cells.
        means.append(lo+s*(hi-lo))
    Q=np.concatenate(Q); W=np.concatenate(W)
    N=np.array(N)
    _,box,widths=envelopes(np.array(f),np.array(E),np.array(means),ref)
    low,high=[fsum(N*z[:,2]) for z in box]
    I=lambda m: fsum(W/np.hypot(Q,m))
    if not low-4e-14*low<=I(ref)<=high+4e-14*high:
        raise ArithmeticError('Reconstructed fixed shape escaped its envelope at reference mass')
    b0=1e-9*x
    lam=d.ds.lam
    if not lam>0:
        raise ArithmeticError('This monotone Hartree envelope requires positive lambda')
    root=lambda m: m*m-b0*b0-lam/2*I(m)
    lo=min(ref,float(np.sqrt(b0*b0+lam*low/2)))
    hi=max(ref,float(np.sqrt(b0*b0+lam*high/2)))
    # A roundoff margin validates the numerical solve; it is not directed rounding.
    margin=5e-12*ref
    full=brentq(root,lo-margin,hi+margin,xtol=max(1e-24,1e-14*ref),rtol=1e-14)
    inverse_exact=fsum(vol/np.hypot(d.ds.p,ref))
    if not lo-margin<=full<=hi+margin or not low-4e-14*low<=inverse_exact<=high+4e-14*high:
        raise ArithmeticError('Nodal reference or continuous representative left the moment root envelope')
    return dict(x=x,grid=g,lifted_histogram=lift,reference_nodal_mass=ref,
                reconstructed_frozen_shape_root=float(full),induced_relative_root_shift=float(full/ref-1),
                admissible_root_range=[lo,hi],admissible_root_relative_width=float((hi-lo)/ref),
                inverse_moment_range=[low,high],inverse_moment_relative_width=fsum(N*widths)/inverse_exact,
                numeric_margin_relative=margin/ref)


def run(a,b):
    return dict(cases=[case(a,False),case(b,True)],
                scope='Frozen positive shapes matching moments at one reference mass; no evolving or physical gap envelope.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('smoke_gauss',type=Path)
    p.add_argument('frozen_scalar',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.smoke_gauss),read(a.frozen_scalar)),indent=2,allow_nan=False))
