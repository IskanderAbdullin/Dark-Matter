"""Common hard thermal mass vs differential shifts for a neutral Yukawa scalar."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_decay_mode_clock_reference import Q, intervals


def massless_normal(p,g):
    # Massless QED normal branch of the HTL pole. This is used only
    # to audit a common hard asymptotic energy shift, not a massive rate.
    f=.125*g*g
    if f==0:return 0.
    hi=f/p
    def h(a):
        return a-hi*(1-a/(2*p)*np.log1p(2*p/a))
    return brentq(h,max(hi*1e-13,1e-40),hi,xtol=1e-16)


def one(s,clock_mass=86.):
    # Set a neutral real scalar coupled off diagonally to two Dirac fields
    # with equal U(1) charge, so only a common thermal-mass shift is used
    # for the gauge branch. Separate shifts below are sensitivity probes.
    rb=s['mass_ratio']
    records=[]
    for r in s['checkpoints']:
        x=r['x'];ms=r['scalar_pole_over_T']
        if x>.101: continue
        vac=intervals(x,ms,rb*x,Q)
        if vac is None: continue
        sv,lo,hi,vw,_=vac
        ev=(hi+lo)/2
        rows=[]
        h=np.sqrt(np.pi**2*106.75/90)*(clock_mass/x)/2.435e18
        for g in (.0003,.001,.003,.01,.03,.10,.30):
            dc=.25*g*g
            tr=g**4*np.log(1/g)/h
            pole=g*g*np.log(1/g)/(4*np.pi)
            for tilt in (0.,.01):
                # Alter the difference while holding the parent's asymptotic contribution.
                P=np.hypot(x,np.sqrt(dc))
                B=np.hypot(rb*x,np.sqrt((1+tilt)*dc))
                e=intervals(P,ms,B,Q)
                init=dict(gauge_coupling=g,daughter_extra_fraction_of_common_mass_sq=tilt,
                          dimensionless_transport_scaling_proxy_over_reference_H=tr,
                          single_line_width_scaling_proxy_over_T=pole)
                if e is None:
                    rows.append(dict(kinematic_status='closed',**init))
                    continue
                E,l,u,span,_=e
                c=.5*(l+u)
                p=np.sqrt(max((c-P)*(c+P),0.))
                omega0=np.hypot(p-Q,B)-np.hypot(p+Q,B)
                assert omega0<0
                displacement=(c-ev)/(vw*.5)
                w=span/vw
                dh=massless_normal(p,g)
                exact=(g*g/4)/(np.hypot(p,g/2)+p)
                rows.append(dict(kinematic_status='open',parent_momentum_center_over_T=p,
                                center_shift_in_vacuum_half_spans=displacement,
                                parent_energy_width_ratio_to_no_thermal_shift=w,
                                single_line_width_proxy_over_daughter_angular_energy_split=pole/abs(omega0),
                                relative_massless_HTL_shift_error_for_asymptotic_dispersion=(exact/dh-1),
                                **init))
        # Shared squared mass delta leaves A unshifted; this independent
        # analytical threshold follows from the exact interval discriminant.
        A=.5*(x*x-(rb*x)**2+ms*ms)
        max_common_sq=(A/ms)**2-x*x
        close=2*np.sqrt(max(max_common_sq,0.))
        records.append(dict(x=x,scalar_pole_over_T=ms,daughter_vacuum_ratio=rb,
                            no_thermal_parent_window=[lo,hi],
                            kinematic_common_gauge_coupling_closure=close,rows=rows))
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    o={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        o[f.name]=one(json.loads(raw))
    print(json.dumps(o,indent=2,allow_nan=False))
