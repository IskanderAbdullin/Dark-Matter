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
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels

def run(grid,duration,count,pair,method="backward",strength=1.,projection="steps",aperture_grid=None,screening_update="frozen",response_smear=.12,compton_strength=0.,same_soft_strength=0.,conversion_strength=0.,snapshots=False,scalar_bins=8,aperture_subdivision=1,witness_mix=0.,response_radial_order=None):
    if duration<=0. or count<1:
        raise ValueError("positive duration and step count required")
    if strength<0. or same_soft_strength<0.:
        raise ValueError("elastic strength cannot be negative")
    if not 0. <= witness_mix < 1. or (witness_mix and projection!="steps"):
        raise ValueError("witness mix is for archived step inputs and must be in [0,1)")
    if scalar_bins<1 or aperture_subdivision<1:
        raise ValueError("positive scalar bins and aperture subdivisions required")
    M=Screened(grid,1,scalar_bins,duration,(1,1,2),.08,strength,aperture_grid,projection,
               subdivision=aperture_subdivision)
    init=M.M.init
    if witness_mix:
        from early_pole_live_aperture_witness_reprojection_reference import initial_projection
        from early_pole_native_aperture_witness_reference import PROFILES
        profiles=json.loads(PROFILES.read_text())
        chosen=next(p for p in profiles if p['grid']==grid)
        _,_,Y0,Y1,*_=initial_projection(chosen,aperture_subdivision,scalar_bins,duration,model=M)
        init=(1-witness_mix)*Y0+witness_mix*Y1
    photon_scattering=None if compton_strength==0. else ComptonKernel(M.M,strength=compton_strength)
    conversion=None if conversion_strength==0. else ConversionKernel(M.M,strength=conversion_strength)
    elastic=M if same_soft_strength==0. else ElasticChannels(M,grid,same_soft_strength)
    active=None if screening_update=="frozen" or not strength else ActiveResponse(M,response_smear,radial_order=response_radial_order)
    updates=[]
    z0=np.r_[init.ravel(),M.ns,0.,0.]
    def total(z):
        _,_,_,r,s=M.M.weights(z[:M.j0].reshape(-1,2),True)
        _,rho,ent,_,_=M.scalar(z[M.j0:M.j0+M.n])
        return r+rho,s+ent
    slow_calls=0
    def advance(z,dt):
        def f(_,y):
            nonlocal slow_calls
            slow_calls+=1
            out=Cold.rhs(M,y,pair)
            if photon_scattering is not None or conversion is not None:
                state=M.M.weights(y[:M.j0].reshape(-1,2),True)[:3]
                for addition in (photon_scattering,conversion):
                    if addition is None:
                        continue
                    A,S,W,*_=addition.collision(*state)
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
    pictures=[]
    def inspect(k,y):
        Y=y[:M.j0].reshape(-1,2)
        e,F,E,R,S=M.M.weights(Y,True)
        pr=M.j0
        n=y[pr:pr+M.n]
        d=Cold.rhs(M,y,False)[pr:pr+M.n]
        pictures.append(dict(step=k,exposure=duration*k/count,cold_bins=n.tolist(),
                             instantaneous_source_by_q=d.tolist(),
                             charged_and_photon_counts=[float(Y[a,0].sum()) for a in M.M.r],
                             hard_logits=F.tolist(),hard_moments=Y.tolist()))
    if snapshots:
        inspect(0,z)
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
        if snapshots and (j+1 in (1,2,4,8,count//4,count//2,3*count//4,count)):
            inspect(j+1,z)
    et,st=total(z)
    assert abs(et-E0)<3e-12 and abs(st-S0-z[-2])<3e-11
    return dict(grid=grid,pair=pair,method=method,strength=strength,projection=projection,
                compton_strength=compton_strength,same_soft_strength=same_soft_strength,
                conversion_strength=conversion_strength,slow_rhs_calls=slow_calls,
                aperture_grid=aperture_grid or grid,duration=duration,steps=count,scalar_gain=states[-1]-states[0],
                state_history=states,energy_drift=et-E0,entropy_closure=st-S0-z[-2],elastic=controls,
                screening_update=screening_update,response_smear=response_smear,response_history=updates,
                maximum_compton_remap_rejection=None if photon_scattering is None else photon_scattering.reject,
                maximum_conversion_remap_rejection=None if conversion is None else conversion.reject,
                hard_edge_sets=M.M.edge_sets if snapshots else None, snapshots=pictures,
                **(dict(scalar_bins=scalar_bins,aperture_subdivision=aperture_subdivision,witness_mix=witness_mix)
                   if scalar_bins!=8 or aperture_subdivision!=1 or witness_mix else {}),
                **(dict(response_radial_order=response_radial_order) if response_radial_order is not None else {}))

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
    p.add_argument("--conversion-strength",type=float,default=0.)
    p.add_argument("--snapshots",action="store_true")
    p.add_argument("--scalar-bins",type=int,default=8)
    p.add_argument("--aperture-subdivision",type=int,default=1)
    p.add_argument("--witness-mix",type=float,default=0.)
    p.add_argument("--response-radial-order",type=int)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
