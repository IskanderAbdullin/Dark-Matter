"""Frozen timelike dispersion only; no occupation of collective modes assumed."""
import json
import numpy as np
from scipy.optimize import brentq

from early_pole_frozen_dynamic_vlasov_reference import Snapshot,load,THERMAL,real_kernels


def eq(x):
    L,T,_,_=real_kernels(x,1.)
    return np.array([L,T])*THERMAL


def compute(fun,plasma):
    ans=[]
    # Only finite-velocity distributions have a finite light-cone value.
    endpoint=None if fun is eq else np.sqrt(-fun(1.)[0])
    for z in (.01,.1,.3,.5,1.,1.5,3.,8.,20.):
        # k normalized to g T_ref; the equilibrium mD_0 is sqrt(2/3) in these units.
        retL=lambda y:z*z+fun(1+np.exp(y))[0]
        retT=lambda y:-z*z*np.exp(y)*(2+np.exp(y))+fun(1+np.exp(y))[1]
        # The thermal longitudinal root redshifts to an exponentially small
        # light-cone gap. For a massive bath it may cease to be timelike.
        eL=retL(-28.)
        # Report only a root that is separated from the floating-point edge.
        xL=1+np.exp(brentq(retL,-28.,12.)) if eL<0. else None
        xT=1+np.exp(brentq(retT,-28.,12.))
        ans.append(dict(momentum=z,omega_transverse=xT*z,
                         transverse_excess_mass_sq=(xT*xT-1)*z*z,
                         omega_longitudinal=xL*z if xL is not None else None,
                         longitudinal_inverse_at_lightcone=eL,
                         k_timelike_endpoint_formula=endpoint))
    assert max(abs(ans[0]['omega_transverse']/np.sqrt(plasma)-1),
               abs(ans[0]['omega_longitudinal']/np.sqrt(plasma)-1))<.002
    return ans


def run():
    rows={'two_thermal':compute(eq,THERMAL/3)}
    for grid in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        for sm in (0.,.18):
            s=Snapshot(d,18,sm)
            rows[grid+str(sm)]=compute(lambda x:s.retarded(x)[:2],s.plasma_frequency_sq())
    return dict(scope='Frozen collisionless retarded dispersion; no imaginary pole widths, '
                      'statistical pole occupation, hard matching or feedback.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
