"""Condition cold coefficients with frozen charges, scalar occupation and residue."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_point_factorial import K,prepare
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from finite_mass_qed_timelike_pole_angular_reference import angle


def estimate(centers,solve):
    rows=[]
    h=2.5e-5
    for k in (0,1):
        def at(s):
            pars=np.array(centers,dtype=float)
            pars[k]*=np.exp(s)
            return solve(*pars)
        raw=[[at(-d),at(d)] for d in (h,h/2)]
        c=np.asarray([np.log(v/u)/(2*d) for (u,v),d in zip(raw,(h,h/2),strict=True)])
        ext=(4*c[1]-c[0])/3
        if np.max(abs(c[1]-c[0])/np.maximum(abs(ext),1e-90))>2e-5 or min(abs(ext))<=0:
            raise ArithmeticError(('Differential collision sensitivity has not resolved',k,c.tolist(),ext.tolist()))
        rows.append(dict(step_sizes_log=(h,h/2),centered_slopes=c.tolist(),richardson_slopes=ext.tolist(),
                         local_log_parameter_change_for_one_percent_coefficient=(np.log1p(.01)/abs(ext)).tolist(),
                         paired_slopes_relative_defect=(abs(c[1]-c[0])/abs(ext)).tolist()))
    return rows


def run(card,thresholds):
    m,(omega,R),(charged,cb),ed,n=prepare(card,2)
    t=[x for x in thresholds if x['grid']=='fine' and x['subcells']==2]
    if len(t)!=1 or t[0]['scalar_mass']!=m:
        raise ArithmeticError('Need the paired numerical threshold enclosure')
    wlo,wref,whi=t[0]['frequencies']
    if abs(wref/omega-1)>2e-12 or not wlo<omega<whi:
        raise ArithmeticError('The saved enclosure does not match the held transverse pole')
    out=[]
    for q in (3.5e-6,4.7e-6):
        f=n[np.searchsorted(ed,q,side='right')-1]
        emax=[]
        def call(mu,w):
            E=np.hypot(q,mu)
            v=[]
            for radial,az in ((10,18),(16,24)):
                a,b,d=angle(charged,K,w,R,E,'T',8,(radial,az),q,
                            occupation_callbacks=cb,mirror_azimuth=True)[0]
                if min(a,b,d)<=0 or abs(a-b-d)>4e-13*max(a,b,d):
                    raise ArithmeticError('Invalid bare point coefficients')
                v.append(q*q/E*np.array([(1+f)*a,f*b]))
            err=np.max(abs(v[0]/v[1]-1))
            if err>1e-7:
                raise ArithmeticError(('Paired point quadratures disagree',q,mu,w,err))
            emax.append(err)
            return v[1]
        rows=estimate((m,omega),call)
        origin=call(m,omega)
        sample=[(call(m,w)/origin-1).tolist() for w in (wlo,whi)]
        out.append(dict(q=q,scalar_occupation=f,base_forward_inverse=origin.tolist(),
                        variations=rows,relative_change_at_frequency_endpoints_holding_residue=sample,
                        maximum_paired_quadrature_defect=max(emax)))
    return dict(x=.05,k=K,axes=('scalar_pole_mass','photon_energy'),columns=('gross_forward','inverse'),
                fixed_residue=R,base_mass=m,base_energy=omega,
                sampled_frequency_interval=[wlo,whi],points=out,
                scope='Local perturbation in either fitted pole parameter at held occupations and residue. The 1% tolerances are local first-order scales, not global error certificates.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('d2',type=Path)
    p.add_argument('thresholds',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.d2),read(a.thresholds)),indent=2,allow_nan=False))
