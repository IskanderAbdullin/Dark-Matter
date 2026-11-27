"""Join the finite hard same-flavor sector to the other slow channels."""
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_same_flavor_hard_redistribution_reference import SameFlavorKernel


def run(grid,duration,level):
    if not duration>0.:
        raise ValueError('duration must be positive')
    C=Cold(grid,1,0,8,True,exposure=duration)
    M=C.M
    order=[(1,3,4,1),(2,4,5,2),(1,4,5,2)][level]
    CH, CV=ComptonKernel(M), ConversionKernel(M)
    K=SameFlavorKernel(M,order,'both',.15)
    z0=np.r_[M.init.ravel(),C.ns,0.,0.]
    a=slice(C.j0,C.j0+C.n)
    def totals(z):
        _,_,_,energy,S=M.weights(z[:C.j0].reshape(-1,2),True)
        _,rho,ent,_,_=C.scalar(z[a])
        return energy+rho,S+ent
    e0,s0=totals(z0)
    t=np.linspace(0.,duration,9)
    result=[]
    for enabled in (False,True):
        def rhs(_,z):
            R=Cold.rhs(C,z,True)
            state=M.weights(z[:C.j0].reshape(-1,2),True)[:3]
            for sector in (CH,CV):
                D,S,W,*_=sector.collision(*state)
                R[:C.j0]+=D.ravel()
                R[-2:]+=[S,W]
            if enabled:
                D,S,W=K.collision(*state)
                R[:C.j0]+=D.ravel()
                R[-2:]+=[sum(S),sum(W)]
            return R
        h=solve_ivp(rhs,(0.,duration),z0,method='DOP853',rtol=2e-10,atol=3e-13,
                    t_eval=t,max_step=duration/4)
        assert h.success
        norms=np.array([totals(z) for z in h.y.T])
        dr=norms[:,0]-e0
        ds=norms[:,1]-s0-h.y[-2]
        assert np.max(abs(dr))<3e-12 and np.max(abs(ds))<3e-11
        gains=(h.y[a]-z0[a,None]).sum(axis=0)
        count_change=[float((h.y[:C.j0,-1]-z0[:C.j0]).reshape(-1,2)[r,0].sum()) for r in M.r]
        assert abs(sum(count_change))<2e-13
        result.append(dict(same_hard=enabled,scalar_gain=gains.tolist(),energy_error=float(np.max(abs(dr))),
                           entropy_error=float(np.max(abs(ds))),rhs_calls=h.nfev,
                           charged_pair_and_photon_count_changes=count_change))
    return dict(grid=grid,aperture_floor=.15,hard_quadrature=order,exposure=t.tolist(),cases=result,
                hard_packets=[g[-1].size for g in K.graph],same_hard_remap_rejection=K.reject)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--duration',type=float,default=2e-6)
    p.add_argument('--level',type=int,choices=(0,1,2),default=0)
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.duration,a.level),indent=2,allow_nan=False))
