"""Conservative split: cold/pair half steps surrounding admissible elastic steps."""
import argparse
import json


import numpy as np
from scipy.integrate import solve_ivp
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_screened_elastic_discrete_gradient_reference import step
from early_pole_screened_active_response_update_reference import ActiveResponse
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels

def run(grid,duration,count,pair,method="backward",strength=1.,projection="steps",aperture_grid=None,screening_update="frozen",response_smear=.12,compton_strength=0.,same_soft_strength=0.):
    if duration<=0. or count<1:
        raise ValueError("positive duration and step count required")
    if strength<0. or same_soft_strength<0.:
        raise ValueError("elastic strength cannot be negative")
    M=Screened(grid,1,8,duration,(1,1,2),.08,strength,aperture_grid,projection)
    photon_scattering=None if compton_strength==0. else ComptonKernel(M.M,strength=compton_strength)
    elastic=M if same_soft_strength==0. else ElasticChannels(M,grid,same_soft_strength)
    active=None if screening_update=="frozen" or not strength else ActiveResponse(M,response_smear)
    updates=[]
    z0=np.r_[M.M.init.ravel(),M.ns,0.,0.]
    def total(z):
        _,_,_,r,s=M.M.weights(z[:M.j0].reshape(-1,2),True)
        _,rho,ent,_,_=M.scalar(z[M.j0:M.j0+M.n])
        return r+rho,s+ent
    def advance(z,dt):
        def f(_,y):
            out=Cold.rhs(M,y,pair)
            if photon_scattering is not None:
                A,S,W,_=photon_scattering.collision(*M.M.weights(y[:M.j0].reshape(-1,2),True)[:3])
                out[:M.j0]+=A.ravel()
                out[-2:]+=[S,W]
            return out
        h=solve_ivp(f,(0.,dt),z,method='DOP853',rtol=2e-10,atol=3e-13)
        assert h.success
        return h.y[:,-1]
    h=duration/count
    z=z0.copy()
    E0,S0=total(z)
    controls=[]
    states=[float(M.ns.sum())]
    for j in range(count):
        z=advance(z,h/2)
        if strength or same_soft_strength:
            Y=z[:M.j0].reshape(-1,2)
            if active is not None:
                updates.append(active.refresh(Y,screening_update))
            f,r=step(elastic,Y,h,method)
            z[:M.j0]=f.ravel()
            z[-2]+=r['entropy_gain']
            controls.append(r)
        z=advance(z,h/2)
        states.append(float(z[M.j0:M.j0+M.n].sum()))
    et,st=total(z)
    assert abs(et-E0)<3e-12 and abs(st-S0-z[-2])<3e-11
    return dict(grid=grid,pair=pair,method=method,strength=strength,projection=projection,
                compton_strength=compton_strength,same_soft_strength=same_soft_strength,
                aperture_grid=aperture_grid or grid,duration=duration,steps=count,scalar_gain=states[-1]-states[0],
                state_history=states,energy_drift=et-E0,entropy_closure=st-S0-z[-2],elastic=controls,
                screening_update=screening_update,response_smear=response_smear,response_history=updates,
                maximum_compton_remap_rejection=None if photon_scattering is None else photon_scattering.reject)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--duration',type=float,default=8e-8)
    p.add_argument('--count',type=int,default=1)
    p.add_argument('--pair',action='store_true')
    p.add_argument("--method",choices=("backward","midpoint"),default="backward")
    p.add_argument("--strength",type=float,default=1.)
    p.add_argument("--projection",choices=("steps","cell-conservative"),default="steps")
    p.add_argument("--aperture-grid",choices=("fine","third"))
    p.add_argument("--screening-update",choices=("frozen","debye","full"),default="frozen")
    p.add_argument("--response-smear",type=float,default=.12)
    p.add_argument("--compton-strength",type=float,default=0.)
    p.add_argument("--same-soft-strength",type=float,default=0.)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
