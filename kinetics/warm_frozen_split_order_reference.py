"""Swapped symmetric pair/collision schedules on the same frozen kernels."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from coupled_warm_pair_continuation_reference import Collision,Late,collide,setup
from warm_stiff_microstep_refinement_reference import sdirk2
from joint_resonant_frw_split_reference import pack,redshift
from late_radiation_background_reference import TabulatedRadiation,MPL
from scalar_pair_expanding_continuation_reference import MASS,moments,rate,reconstruct,e_total

def run(eos,seed,h,methods,schemes,steps,pair_strength=1.):
    assert 0.<pair_strength<=1.
    s=json.loads(Path(seed).read_text())["snapshots"][-1]
    q=np.asarray(s["lab_momentum_edges_GeV"])
    th=[np.asarray(s[x]) for x in ("entropy_variable_scalar","entropy_variable_fermion")]
    b=np.exp(h/2)
    p=q/b
    initial,*_=redshift(q,1.,b,th)
    T0=s["T_visible_GeV"]
    clock=TabulatedRadiation(eos,T0)
    tvis,rhov,*_=clock.state(b)
    U=moments(p,initial); E=e_total(U)
    H=np.sqrt((rhov+E)/3.)/MPL
    k=Collision(p)
    assert len(k.y)
    pilot=Late(k,H,[],True)
    u=pack(initial)
    _,_,T=pilot.reference(u,tvis)
    panels=[]
    setup(max(T,tvis),k,T,pilot.eta,3,2,2,8,32.,.82,cache=panels,cache_only=True)
    model=Late(k,H,panels,True,compact=True,local_jac=True)
    N=U[:,0].sum(axis=1)
    def evaluate(order,n):
        x=u.copy()
        total_out=np.zeros((2,2))
        fail=0;accepted=0
        def pair(y,dt):
            dt*=pair_strength
            mom=moments(p,model.arrays(y))
            F,_,_=rate(p,model.arrays(y),(3,7,6,2))
            mid=reconstruct(p,mom+(dt/(2*H))*F[:,:,1:-1])
            F,_,_=rate(p,mid,(3,7,6,2))
            V=F*dt/H
            total_out[:]+=V[:,:,[0,-1]].sum(axis=2)
            return pack(reconstruct(p,mom+V[:,:,1:-1]))
        def con(y,t):
            nonlocal fail,accepted
            if "sdirk2" in methods:
                z,_=sdirk2(model,y,t,T); acc,rej=1,0
            else:
                z,r,acc,rej=collide(model,y,t,T)
            fail+=rej
            accepted+=acc
            return z
        for _ in range(n):
            d=h/n
            if order=="pair-collision-pair":
                x=pair(con(pair(x,d/2),d),d/2)
            elif order=="pair-two-half-collisions-pair":
                x=pair(con(con(pair(x,d/2),d/2),d/2),d/2)
            else:
                x=con(pair(con(x,d/2),d),d/2)
        end=moments(p,model.arrays(x)).sum(axis=2)+total_out
        c=end[:,0]
        e=float(end[:,1].sum()+MASS@c)
        assert abs(e/E-1.)<2e-9 and abs(c.sum()/N.sum()-1.)<2e-9
        return dict(outer_subintervals=n,count_relative_to_initial=(c/N).tolist(),energy_ledger=e/E-1.,
                    count_ledger=c.sum()/N.sum()-1.,collision_accept_reject=[accepted,fail])
    out={}
    for key in schemes:
        out[key]=[]
        for n in steps:
            out[key].append(evaluate(key,n))
            print(key,n,file=sys.stderr,flush=True)
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("seed");p.add_argument("--step",type=float,default=.006)
    p.add_argument("--method",choices=("sdirk2","backward"),default="backward")
    scheme=("pair-collision-pair","collision-pair-collision","pair-two-half-collisions-pair")
    p.add_argument("--schemes",nargs="+",choices=scheme,default=scheme)
    p.add_argument("--substeps",nargs="+",type=int,default=[1,2,4])
    p.add_argument("--pair-strength",type=float,default=1.)
    a=p.parse_args();print(json.dumps(run(a.eos,a.seed,a.step,a.method,a.schemes,a.substeps,a.pair_strength),indent=2,allow_nan=False))
