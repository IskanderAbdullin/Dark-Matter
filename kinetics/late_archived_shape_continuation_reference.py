"""Source-free continuation straight from the archived 12 MeV entropy shells."""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import critical_T,from_row,fugacity,moments,solve
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_combined_centered_implicit_reference import Combined
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from late_soft_two_moment_test_reference import totals
from late_two_temperature_expansion_reference import inverse
from three_species_decay_onsager_reference import m,mp,P
from sm_dark_entropy_exchange_reference import StandardModel,MPL


class Late(Combined):
    def __init__(self,k,H,panels,prompt):
        super().__init__(k,H,panels)
        self.prompt=prompt

    def reference(self,q,guess):
        U,_,_=self.moment(q)
        V=U.reshape(2,self.n,2)
        n=V[:,:,0].sum(axis=1)
        self.t=max(guess,1.01*critical_T(n[0],20))
        if self.prompt:
            K=sum((V[j,:,0]*c.low+V[j,:,1]*c.dk).sum() for j,c in enumerate(self.k.c))
            normal=solve(n,K)
            if normal["normal_common_phase"]: self.t=normal["common_T_GeV"]
        et=[fugacity(self.t,c.mass,n[j],j==0,20) for j,c in enumerate(self.k.c)]
        self.eta=[c.mass/self.t+et[j] for j,c in enumerate(self.k.c)]
        out=pack([np.array([c.low/self.t-et[j],c.dk/self.t]) for j,c in enumerate(self.k.c)])
        R0=np.exp(sum(self.eta)-self.k.E[self.k.y]/self.t) if self.prompt else None
        return out,R0,self.t

    def centered(self,z,ref,R0):
        return super().centered(z,ref,R0) if self.prompt else self.elastic(z)


def collide(model,u,h,T,depth=0):
    errors=[]
    empty=np.any(model.moment(u)[0][::2]==0.)
    guesses=(True,) if empty else ((True,False) if model.prompt else (False,True))
    for ref in guesses:
        try:
            v,l=advance(model,u,h,T,reference_guess=ref,vacuum_scaled=True)
            return v,l["residual"],1,len(errors)
        except (RuntimeError,FloatingPointError,RuntimeWarning) as e:
            errors.append(e)
    if depth==5: raise errors[-1]
    v,r1,n1,f1=collide(model,u,h/2,T,depth+1)
    w,r2,n2,f2=collide(model,v,h/2,T,depth+1)
    return w,max(r1,r2),n1+n2,f1+f2+len(errors)


def report(model,th,a,Ninit):
    U,S,_=model.moment(pack(th))
    N,E=totals(model.k,U,a)
    K=(E/a-N*np.array([mp,m]))/a**3
    n=N/a**3
    Ts=[inverse(x,y,j,"K",.010/a)["T"] for j,(x,y) in enumerate(zip(n,K))]
    scalar_fraction=n[0]/moments(Ts[0],mp,0.,True)[0]
    v=U.reshape(2,model.n,2)
    theta_min=np.min(th[0][0]+np.minimum(0.,th[0][1]))
    assert theta_min>0.
    maximum= np.exp(-theta_min)/(-np.expm1(-theta_min))
    return dict(b=a,equivalent_T_GeV=Ts,scalar_count_over_own_saturation=scalar_fraction,
                maximum_resolved_scalar_occupation=float(maximum),
                scalar_first_cell_number_fraction=v[0,0,0]/v[0,:,0].sum(),
                count_defect=float(max(abs(N/Ninit-1.))), dark_entropy_a3=S*a**3),E


def subdivide(q,th,stride):
    """Exact subdivision of a piecewise-linear kinetic entropy variable."""
    if stride<=0:return q,th
    sel=np.arange(0,len(q)-1,stride)
    r=np.sort(np.r_[q,(q[sel]+q[sel+1])/2])
    idx=np.searchsorted(q,(r[:-1]+r[1:])/2,"right")-1
    out=[]
    for mass,v in zip((mp,m),th):
        k=lambda p:p*p/(np.hypot(p,mass)+mass)
        a,b=k(q[:-1])[idx],np.diff(k(q))[idx]
        out.append(np.stack([v[0,idx]+v[1,idx]*(k(r[:-1])-a)/b,
                             v[1,idx]*np.diff(k(r))/b]))
    return r,out

def run(eos,row,snap,nx=12,stop=410.,prompt_until=8.,orders=(3,2,2,8),off=False,radiation_table=None,dark_gravity=False,refinement=0,preroll_stop=None,preroll_steps=0):
    assert not off or prompt_until==0.
    assert not dark_gravity or radiation_table is not None
    assert row["nsteps"]==snap["steps"]
    base=from_row(row,eos,20)
    tab=None
    if radiation_table is not None:
        from late_radiation_background_reference import TabulatedRadiation
        tab=TabulatedRadiation(radiation_table,base["actual_visible_T_GeV"])
    def H(a,rhod=0.):
        y=(np.pi**2/15)*(base["actual_visible_T_GeV"]/a)**4 if tab is None else tab.state(a)[1]
        return np.sqrt((y+(rhod if dark_gravity else 0.))/3)/MPL
    q=np.asarray(snap["lab_momentum_edges_GeV"])
    th=[np.asarray(snap["entropy_variable_"+j]) for j in ("scalar","fermion")]
    initial_first_boundary=q[1]
    original_q=q.copy()
    q,th=subdivide(q,th,refinement)
    ks=Collision(q); model=Late(ks,H(1.),[],False)
    N0,E0=totals(ks,model.moment(pack(th))[0],1.)
    assert max(abs(N0/np.asarray(base["species_counts_GeV3"])-1.))<4e-12
    init,_=report(model,th,1.,N0)
    # Use identical cold cutoff after subdivision; sum whole exact bins.
    def cold_fraction(mod,t):
        U=mod.moment(pack(t))[0].reshape(2,mod.n,2)
        return float(U[0,q[1:]<=initial_first_boundary*(1+1e-12),0].sum()/U[0,:,0].sum())
    init["scalar_count_below_archived_first_edge"]=cold_fraction(model,th)
    e0=E0.sum()-np.array([mp,m])@N0
    assert abs(e0/base["retained_dark_kinetic_energy_GeV4"]-1.)<4e-11
    pressure=np.zeros(2); Q4=np.zeros(2); Q3=np.zeros(2); maxres=0.; accepted=0; rejected=0
    steps=[]
    prompt_geometric_end=q[-1]/P
    cut=min(prompt_until,prompt_geometric_end)
    if preroll_stop is not None:
        assert 1.<preroll_stop<stop and preroll_steps>0
        cut0=np.log(preroll_stop)
        grid=np.r_[np.linspace(0.,cut0,preroll_steps+1),np.linspace(cut0,np.log(stop),nx+1)[1:]]
    else:
        grid=np.linspace(0.,np.log(stop),nx+1)
    if 1.<cut<stop: grid=np.unique(np.r_[grid,np.log(cut)])
    n_prompt=0
    for x,y in zip(grid[:-1],grid[1:]):
        a,b,c=np.exp([x,(x+y)/2,y])
        th,t,*_=redshift(q,a,b,th); pressure+=t
        if not off:
            kk=Collision(q/b)
            use_prompt=(b<cut and len(kk.y)>0)
            n_prompt+=int(use_prompt)
            if b>=prompt_geometric_end*(1+1e-9):
                assert len(kk.y)==0
            local=Late(kk,H(b),[],False)
            u=pack(th)
            E_mid=totals(kk,local.moment(u)[0],b)[1]
            hb=H(b,E_mid.sum()/b**4)
            first=Late(kk,hb,[],use_prompt)
            _,_,T=first.reference(u,base["common_T_GeV"]/b)
            panels=[]
            setup(max(T,base["common_T_GeV"]/b),kk,T,first.eta,*orders,32.,.82,cache=panels)
            co=Late(kk,hb,panels,use_prompt)
            uu,r,acc,rej=collide(co,u,y-x,T)
            maxres=max(maxres,r);accepted+=acc;rejected+=rej
            dE=totals(kk,co.moment(uu)[0],b)[1]-totals(kk,co.moment(u)[0],b)[1]
            Q4+=dE;Q3+=dE/b
            th=co.arrays(uu)
        th,t,*_=redshift(q,b,c,th);pressure+=t
        k=Collision(q/c)
        diagnostic=Late(k,H(c),[],False)
        summary,E=report(diagnostic,th,c,N0)
        summary["scalar_count_below_archived_first_edge"]=cold_fraction(diagnostic,th)
        summary.update(H_GeV=float(H(c,E.sum()/c**4)),
                       photon_temperature_GeV=tab.state(c)[0] if tab else base["actual_visible_T_GeV"]/c)
        steps.append(summary)
    ledger=max(abs(E-E0-pressure-Q4))/E0.sum()
    balance=abs(Q4.sum())/E0.sum()
    assert max(ledger,balance,max(x["count_defect"] for x in steps))<5e-9
    assert np.min(np.diff([init["dark_entropy_a3"]]+[x["dark_entropy_a3"] for x in steps]))>=-5e-11*init["dark_entropy_a3"]
    return dict(nsteps=nx,elastic_orders=orders,elastic=not off,prompt_until=prompt_until,
                initial=init,steps=steps,
                preroll_stop=preroll_stop,preroll_steps=preroll_steps,
                upper_support_kinematic_prompt_end=float(prompt_geometric_end),
                intervals_with_prompt=n_prompt,actual_intervals=len(steps),
                background="tabulated radiation" if tab else "freely redshifted photons only",
                includes_retained_dark_energy_in_H=dark_gravity,initial_first_cell_upper_q_over_TR=float(initial_first_boundary*row["final_a"]/row["initial_visible_T_GeV"]),
                initial_archive_shells=len(original_q)-1,shells=len(q)-1,refinement_stride=refinement,
                accepted_collision_solves=accepted,rejected_collision_attempts=rejected,max_collision_residual=maxres,
                energy_ledger_defect=float(ledger),opposite_heat_defect=float(balance),
                pressure_trace_a4_energy=pressure.tolist(),
                collision_energy_a4=Q4.tolist(),
                energy_ledger_defect_in_initial_kinetic_units=float(max(abs(E-E0-pressure-Q4))/e0),
                collision_scalar_a3_heat_over_initial_a3_kinetic=float(Q3[0]/e0))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("run");p.add_argument("shape")
    p.add_argument("--steps",type=int,default=24);p.add_argument("--stop",type=float,default=410.)
    p.add_argument("--prompt-until",type=float,default=8.);p.add_argument("--off",action="store_true")
    p.add_argument("--orders",type=int,nargs=4,default=(3,2,2,8))
    p.add_argument("--radiation-table")
    p.add_argument("--dark-gravity",action="store_true")
    p.add_argument("--refinement-stride",type=int,default=0)
    p.add_argument("--preroll-stop",type=float)
    p.add_argument("--preroll-steps",type=int,default=0)
    a=p.parse_args()
    b=json.loads(Path(a.run).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    print(json.dumps(run(StandardModel(a.eos),b,s,a.steps,a.stop,a.prompt_until,a.orders,a.off,a.radiation_table,a.dark_gravity,a.refinement_stride,a.preroll_stop,a.preroll_steps),indent=2,allow_nan=False))
