"""Independently re-evaluate CLASS momentum-grid shape convergence."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp


def sample(x,f,qcap,n,transform=False):
    P=CubicSpline(x,f)
    if transform:
        t=np.arange(1,n+1)/(n+1.)
        q=1/t-1
        tail=np.maximum(0.,q-x[-1])*(f[-1]-f[-2])/(f[-1]*(x[-1]-x[-2]))
        w=P(np.clip(q,x[0],x[-1]))*np.exp(tail)/(t*t*(n+1))
        spacing=None
    else:
        q=np.linspace(qcap/n,qcap,n)
        w=np.where(q<x[0],f[0],P(q))*(qcap/n)
        assert q[-1]<=x[-1]
        w[-1]*=.5
        spacing=qcap/n
    assert np.min(w)>0.
    iq=np.array([np.sum(w*q**k) for k in (2,3,4)])
    return iq/iq[0],spacing



def run(eos,paths):
    b=Background(eos)
    Wgam=b.ogamma*(b.H0/100)**2
    coeff=b.T0/(mp*Wgam*(7/8)*(4/11)**(4/3))
    d={}
    for root in paths:
        for name in ("direct50_absolute","direct100_absolute"):
            lines=(root/(name+".ini")).read_text().splitlines()
            c=dict((a.strip(),b.strip()) for a,b in (s.split("=",1) for s in lines if "=" in s))
            assert c["ncdm_quadrature_strategy"] in ("2","3")
            transformed=c["ncdm_quadrature_strategy"]=="2"
            x,f=np.loadtxt(c["ncdm_psd_filenames"]).T
            z=np.array([np.trapezoid(f*x**k,x) for k in (2,3,4)])
            moment,spacing=sample(x,f,float(c["ncdm_maximum_q"]),int(c["ncdm_N_momentum_bins"]),transformed)
            out=dict(momentum_grid="1/t-1" if transformed else "linear",momentum_grid_spacing=spacing,mean_q_and_q2_full_table=(z[1:]/z[0]).tolist(),
                     mean_q_and_q2_native_grid=moment[1:].tolist(),
                     fractional_errors_of_mean_q_and_q2=(moment[1:]/(z[1:]/z[0])-1).tolist(),
                     tail_q2_fraction_above_grid_cap=None if transformed else float(np.trapezoid((f*x**4)[x>float(c["ncdm_maximum_q"])],
                                                                            x[x>float(c["ncdm_maximum_q"])])/z[-1]))
            log=root/(name+".log")
            if log.exists() and log.stat().st_size>0:
                native=float(re.search(r"gives Delta N_eff = ([^\s]+)",log.read_text())[1])
                own=coeff*float(c["omega_ncdm"])*moment[1]
                assert abs(own/native-1.)<2e-5
                out.update(class_relativistic_delta_Neff=native,calculated_relativistic_delta_Neff=own)
            postfix="transformed" if transformed else f'cap{c["ncdm_maximum_q"]}'
            d.setdefault(f'{c["ncdm_N_momentum_bins"]}q-{postfix}',{})[name]=out
    return dict(scope="Cubic PSD plus linear or transformed q trapezoid; comparison to continuous input moments "
                      "and the CLASS logged Delta Neff where complete. Neither a numerical guarantee "
                      "on the hierarchy nor a BBN or CMB constraint.", tiers=d)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("folders",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.folders),indent=2,allow_nan=False))
