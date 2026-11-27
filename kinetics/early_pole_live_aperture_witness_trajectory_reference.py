"""Conservative 1↔2 trajectory starting from an explicitly projected witness."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES


def solve(C,Y0,Y1,alpha,pair,refine):
    y0=np.r_[((1-alpha)*Y0+alpha*Y1).ravel(),C.ns,0.,0.]
    r=C.M.r
    s=slice(C.j0,C.j0+C.n)
    def totals(y):
        hard=y[:C.j0].reshape(-1,2)
        _,_,_,e,S=C.M.weights(hard,True)
        _,en,ent,_,m=C.scalar(y[s])
        return np.r_[e+en,S+ent,[hard[z,0].sum() for z in r],y[s].sum(),m]
    exposure=np.r_[0.,np.geomspace(C.exposure/1024,C.exposure,17)]
    def rhs(_,y):return C.rhs(y,pair)
    rr=solve_ivp(rhs,(0.,C.exposure),y0,method='DOP853',rtol=2e-10/refine,
                 atol=3e-13/refine,dense_output=False,t_eval=exposure,max_step=C.exposure/(4*refine))
    if not rr.success:raise RuntimeError(rr.message)
    T=np.array([totals(y) for y in rr.y.T])
    K=(T-T[0]).T
    rates=np.array([C.rhs(y,False)[s].sum() for y in rr.y.T])
    assert np.max(np.abs(K[0]))<3e-12
    assert np.max(np.abs(K[1]-rr.y[-2]))<3e-11
    assert np.min(rr.y[-2])>=-1e-13
    N=K[5]
    if not pair:
        assert max(np.max(abs(K[3]+N)),np.max(abs(K[4]-N)))<5e-13
        assert np.max(abs(K[2]))<5e-13
    else:
        # 2 photons <-> 2 charged particles conserves combined number; decay
        # conserves charged count and adds/removes one scalar independently.
        assert np.max(abs(K[2]+K[3]+K[4]))<5e-13
    return dict(witness_mix=alpha,pair=pair,input_total_energy=T[0,0],
                duration=C.exposure,sample_exposures=exposure.tolist(),cold_gain=N.tolist(),
                instantaneous_cold_source=rates.tolist(),poles=T[:,6].tolist(),
                charge_and_photon_initial_counts=T[0,2:5].tolist(),
                charge_and_photon_changes=K[2:5].T.tolist(),
                max_energy_defect=float(np.max(abs(K[0]))),
                max_entropy_ledger_error=float(np.max(abs(K[1]-rr.y[-2]))),
                rhs_calls=rr.nfev,adaptivity_refinement=refine)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--subdivision',type=int,default=16)
    p.add_argument('--scalar-bins',type=int,default=16)
    p.add_argument('--duration',type=float,default=2e-6)
    p.add_argument('--mix',type=float,nargs='+',default=[0.,.85,.97])
    p.add_argument('--pair',action='store_true')
    p.add_argument('--refine',type=int,default=1)
    a=p.parse_args()
    profiles=json.loads(PROFILES.read_text())
    profile=next(x for x in profiles if x['grid']==a.grid)
    C,_,Y0,Y1,*_=initial_projection(profile,a.subdivision,a.scalar_bins,a.duration)
    out=dict(grid=a.grid,subdivision=a.subdivision,scalar_bins=a.scalar_bins,
             cases=[solve(C,Y0,Y1,x,a.pair,a.refine) for x in a.mix])
    print(json.dumps(out,indent=2,allow_nan=False))
