"""Collision-only driver for the first variation in pair strength."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson,trapezoid
from coupled_warm_pair_continuation_reference import Collision,Late,collide,setup
from joint_resonant_frw_split_reference import pack,redshift
from late_radiation_background_reference import TabulatedRadiation,MPL
from scalar_pair_expanding_continuation_reference import moments,rate,e_total
from warm_stiff_microstep_refinement_reference import sdirk2

GRID=np.r_[0.,np.geomspace(1e-5,.1,13),np.linspace(.1,1.,9)[1:]]

def build(eos,seed,h):
    s=json.loads(Path(seed).read_text())["snapshots"][-1]
    q=np.asarray(s["lab_momentum_edges_GeV"])
    th=[np.asarray(s[x]) for x in ("entropy_variable_scalar","entropy_variable_fermion")]
    b=np.exp(h/2)
    p=q/b
    initial,*_=redshift(q,1.,b,th)
    T0=s["T_visible_GeV"]
    clock=TabulatedRadiation(eos,T0)
    tvis,rhov,*_=clock.state(b)
    U=moments(p,initial)
    H=np.sqrt((rhov+e_total(U))/3.)/MPL
    k=Collision(p)
    pilot=Late(k,H,[],True)
    u=pack(initial)
    _,_,T=pilot.reference(u,tvis)
    panels=[]
    setup(max(T,tvis),k,T,pilot.eta,3,2,2,8,32.,.82,cache=panels,cache_only=True)
    return p,u,U,H,T,Late(k,H,panels,True,compact=True,local_jac=True)

def run(eos,seed,h,method):
    p,u,U,H,T,model=build(eos,seed,h)
    N=U[:,0].sum(axis=1)
    x,y=GRID[:-1],GRID[1:]
    result=[]
    for split in (1,2):
        grid=GRID if split==1 else np.sort(np.r_[GRID,(x+y)/2])
        q=u.copy()
        rates=[]; ri=[]; countmax=0.; energymax=0.
        for i,b in enumerate(grid):
            if i:
                dt=h*(b-grid[i-1])
                if method=="sd":
                    q,_=sdirk2(model,q,dt,T)
                else:
                    q,_,_,_=collide(model,q,dt,T)
            th=model.arrays(q)
            K=moments(p,th)
            d=rate(p,th,(3,7,6,2))
            rates.append((d[0][:,0].sum(axis=1)/H/N).tolist())
            ri.append(d[2])
            countmax=max(countmax,float(max(abs(K[:,0].sum(axis=1)/N-1.))))
            energymax=max(energymax,abs(e_total(K)/e_total(U)-1.))
        print(method,split,len(grid),flush=True,file=sys.stderr)
        a=np.asarray(rates)
        result.append(dict(nodes=len(grid),time_fraction=grid.tolist(),net_rate_per_initial_species_per_log_b=rates,inverse_over_forward=ri,
                     composite_trapezoid_count_delta=trapezoid(a,x=grid,axis=0).tolist(),
                     composite_simpson_count_delta=simpson(a,x=grid,axis=0).tolist(),
                     collision_only_count_and_energy_error=[countmax,energymax]))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("seed")
    p.add_argument("--interval",type=float,default=.006)
    p.add_argument("--method",choices=("sd","be"),default="sd")
    a=p.parse_args()
    r=run(a.eos,a.seed,a.interval,a.method)
    for u in r:
        for k in ("composite_trapezoid_count_delta","composite_simpson_count_delta"):
            u[k]=[x*a.interval for x in u[k]]
    print(json.dumps(r,indent=2,allow_nan=False))
