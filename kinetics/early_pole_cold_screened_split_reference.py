"""Conservative split: cold/pair half steps surrounding admissible elastic steps."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile


import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import solve_ivp
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_screened_elastic_discrete_gradient_reference import step
from early_pole_screened_active_response_update_reference import ActiveResponse
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels

def state_digest(x):
    return hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest()


def atomic_write(path, value):
    path=Path(path)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,prefix=path.name+'.',
                                         suffix='.tmp',delete=False) as fp:
            temporary=Path(fp.name)
            json.dump(value,fp,allow_nan=False,separators=(',',':'))
            fp.flush()
            os.fsync(fp.fileno())
        os.replace(temporary,path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def run(grid,duration,count,pair,method="backward",strength=1.,projection="steps",aperture_grid=None,screening_update="frozen",response_smear=.12,compton_strength=0.,same_soft_strength=0.,conversion_strength=0.,snapshots=False,scalar_bins=8,aperture_subdivision=1,witness_mix=0.,response_radial_order=None,
        checkpoint=None,resume=False,checkpoint_every=1,stop_after=None,progress_every=0,elastic_energy_budget=None,macro_step_weights=None,outer_order='slow-elastic-slow'):
    if elastic_energy_budget is not None and (not np.isfinite(elastic_energy_budget) or elastic_energy_budget<=0. or method!="backward"):
        raise ValueError("A finite positive explicit elastic energy budget requires backward steps.")
    if resume and not checkpoint:
        raise ValueError("resume requires a checkpoint path")
    if checkpoint_every<1 or progress_every<0:
        raise ValueError("invalid checkpoint or progress interval")
    if stop_after is not None and (not checkpoint or not 1<=stop_after<=count):
        raise ValueError("stop-after requires a checkpoint and a valid step")
    if checkpoint and not resume and Path(checkpoint).exists():
        raise FileExistsError("checkpoint exists: use explicit resume or choose a new path")
    if duration<=0. or count<1:
        raise ValueError("positive duration and step count required")
    if strength<0. or same_soft_strength<0.:
        raise ValueError("elastic strength cannot be negative")
    if not 0. <= witness_mix < 1. or (witness_mix and projection!="steps"):
        raise ValueError("witness mix is for archived step inputs and must be in [0,1)")
    if scalar_bins<1 or aperture_subdivision<1:
        raise ValueError("positive scalar bins and aperture subdivisions required")
    if outer_order not in ('slow-elastic-slow','elastic-slow-elastic'):
        raise ValueError('Unsupported outer ordering.')
    if outer_order=='elastic-slow-elastic' and (screening_update!='frozen' or not (strength or same_soft_strength)):
        raise ValueError('The alternate order is tested only with active elastic terms and held response.')
    weights=None
    if macro_step_weights is not None:
        weights=tuple(macro_step_weights)
        if len(weights)!=count or any(not isinstance(x,int) or isinstance(x,bool) or x<1 for x in weights):
            raise ValueError('Expected one positive integer macro weight per interval.')
        unit=duration/sum(weights)
        boundaries=np.r_[0,np.cumsum(weights)]*unit
        boundaries[-1]=duration
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
    if progress_every:
        print("initialized fixed supports and channel graphs",file=sys.stderr,flush=True)
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
        pictures.append(dict(step=k,exposure=duration*k/count if weights is None else float(boundaries[k]),cold_bins=n.tolist(),
                             instantaneous_source_by_q=d.tolist(),
                             charged_and_photon_counts=[float(Y[a,0].sum()) for a in M.M.r],
                             hard_logits=F.tolist(),hard_moments=Y.tolist()))
    physics=dict(grid=grid,duration=duration,count=count,pair=pair,method=method,
                 strength=strength,projection=projection,aperture_grid=aperture_grid,
                 screening_update=screening_update,response_smear=response_smear,
                 compton_strength=compton_strength,same_soft_strength=same_soft_strength,
                 conversion_strength=conversion_strength,snapshots=snapshots,scalar_bins=scalar_bins,
                 aperture_subdivision=aperture_subdivision,witness_mix=witness_mix,
                 response_radial_order=response_radial_order)
    if elastic_energy_budget is not None:
        physics["elastic_energy_budget"]=elastic_energy_budget
    if weights is not None:
        physics['macro_step_weights']=list(weights)
    if outer_order!='slow-elastic-slow':
        physics['outer_order']=outer_order
    if checkpoint:
        physics.update(format=1,numpy=np.__version__,scipy=scipy_version,state_initial=state_digest(z0),
                       graph_code=[hashlib.sha256(Path(__file__).with_name(n+'.py').read_bytes()).hexdigest()
                                   for n in (Path(__file__).stem,'early_pole_archived_cold_screened_elastic_reference',
                                             'early_pole_same_flavor_screened_pole_reference',
                                             'early_pole_screened_active_response_update_reference',
                                             'early_pole_screened_elastic_discrete_gradient_reference',
                                             'early_pole_archived_shared_cold_window_reference',
                                             'early_pole_archived_live_compton_reference',
                                             'early_pole_archived_live_flavor_conversion_reference')])
    start=0
    states=[float(M.ns.sum())]
    if resume:
        old=json.loads(Path(checkpoint).read_text())
        if old.get('physics')!=physics:
            raise ValueError("checkpoint does not match the inputs or solver version")
        z=np.asarray(old['state'],dtype=float)
        start=old['steps_done']
        states=old['states']
        controls=old['controls']
        pictures=old['pictures']
        updates=old['updates']
        slow_calls=old['slow_calls']
        if (z.shape!=z0.shape or not np.all(np.isfinite(z)) or old['state_digest']!=state_digest(z)
                or not isinstance(start,int) or not 0<=start<=count or len(states)!=start+1
                or len(controls)!=(start if strength or same_soft_strength else 0)
                or len(updates)!=(start if active is not None else 0)
                or states[0]!=float(M.ns.sum()) or states[-1]!=float(z[M.j0:M.j0+M.n].sum())):
            raise ValueError("invalid checkpoint state or history")
        if photon_scattering is not None:
            photon_scattering.reject=list(old['compton_reject'])
        if conversion is not None:
            conversion.reject=float(old['conversion_reject'])
    elif snapshots:
        inspect(0,z)
    limit=count if stop_after is None else stop_after
    if limit<start:
        raise ValueError("stop-after precedes the restored step")
    def store(k):
        atomic_write(checkpoint,dict(physics=physics,steps_done=k,state=z.tolist(),
                                     state_digest=state_digest(z),states=states,controls=controls,pictures=pictures,
                                     updates=updates,slow_calls=slow_calls,
                                     compton_reject=None if photon_scattering is None else photon_scattering.reject,
                                     conversion_reject=None if conversion is None else conversion.reject))
    for j in range(start,limit):
        if weights is not None:
            h=weights[j]*unit
        if outer_order=='elastic-slow-elastic':
            def elastic_half(y):
                Y=y[:M.j0].reshape(-1,2)
                f,r=step(elastic,Y,h/2,method,energy_budget=elastic_energy_budget)
                y[:M.j0]=f.ravel()
                y[-2]+=r['entropy_gain']
                return r
            before=elastic_half(z)
            z=advance(z,h)
            after=elastic_half(z)
            parts=[v for q in (before,after) for v in q.get('elastic_substeps',[q])]
            controls.append(dict(h=h,rho_drift=before['rho_drift']+after['rho_drift'],
                                 entropy_gain=before['entropy_gain']+after['entropy_gain'],
                                 elastic_substeps=parts))
        else:
            z=advance(z,h/2)
            if strength or same_soft_strength:
                Y=z[:M.j0].reshape(-1,2)
                if active is not None:
                    updates.append(active.refresh(Y,screening_update))
                f,r=step(elastic,Y,h,method,energy_budget=elastic_energy_budget)
                z[:M.j0]=f.ravel()
                z[-2]+=r['entropy_gain']
                controls.append(r)
            z=advance(z,h/2)
        states.append(float(z[M.j0:M.j0+M.n].sum()))
        if snapshots and (j+1 in (1,2,4,8,count//4,count//2,3*count//4,count)):
            inspect(j+1,z)
        if progress_every and ((j+1)%progress_every==0 or j+1==limit):
            print(f"steps {j+1}/{count}; cold gain {states[-1]-states[0]:.8e}",file=sys.stderr,flush=True)
        if checkpoint and ((j+1)%checkpoint_every==0 or j+1==limit):
            store(j+1)
    if limit<count:
        return dict(checkpoint_only=True,steps_done=limit,total_steps=count,scalar_gain_so_far=states[-1]-states[0])
    et,st=total(z)
    assert abs(et-E0)<3e-12 and abs(st-S0-z[-2])<3e-11
    return dict(grid=grid,pair=pair,method=method,strength=strength,projection=projection,**(dict(macro_step_weights=list(weights),macro_exposure=boundaries.tolist()) if weights is not None else {}),
                **(dict(outer_order=outer_order) if outer_order!='slow-elastic-slow' else {}),
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
                **(dict(response_radial_order=response_radial_order) if response_radial_order is not None else {}),
                **(dict(elastic_energy_budget=elastic_energy_budget) if elastic_energy_budget is not None else {}))

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
    p.add_argument("--elastic-energy-budget",type=float)
    p.add_argument("--checkpoint")
    p.add_argument("--resume",action="store_true")
    p.add_argument("--checkpoint-every",type=int,default=1)
    p.add_argument("--stop-after",type=int)
    p.add_argument("--progress-every",type=int,default=0)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
