"""Source-free continuation straight from the archived 12 MeV entropy shells."""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import critical_T,from_row,fugacity,moments,solve
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_combined_centered_implicit_reference import Combined
from joint_resonant_centered_stiff_reference import Centered
from compact_elastic_panel_reference import compile_panels,compact_elastic
from pivoted_prompt_reference import pivoted
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from late_soft_two_moment_test_reference import totals
from late_two_temperature_expansion_reference import inverse
from three_species_decay_onsager_reference import m,mp,P
from sm_dark_entropy_exchange_reference import StandardModel,MPL


class Late(Combined):
    def __init__(self,k,H,panels,prompt,compact=False,pivot=False):
        super().__init__(k,H,panels)
        self.prompt=prompt
        self.compiled=compile_panels(panels,self.n) if compact else None
        self.pivot=pivot

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

    def elastic(self,z):
        return (super().elastic(z) if self.compiled is None else
                compact_elastic(z,self.t,self.eta,self.H,self.compiled,self.d,True))

    def centered(self,z,ref,R0):
        if not self.prompt: return self.elastic(z)
        a,da=pivoted(self,z,ref,R0) if self.pivot else Centered.centered(self,z,ref,R0)
        if not self.panels: return a,da
        b,db=self.elastic(z)
        return a+b,lambda v:da(v)+db(v)


def collide(model,u,h,T,depth=0,continuation=False,verbose=False):
    def obs(i,r,j):
        if i in (0,1,2,5,10,20,40,80):
            import sys
            print(f"collision h={h:.6g}, iter={i}, resid={r:.6g}, coordinate={j}",file=sys.stderr,flush=True)
    observe=obs if verbose else None
    errors=[]
    empty=np.any(model.moment(u)[0][::2]==0.)
    guesses=(True,) if empty else ((True,False) if model.prompt else (False,True))
    for ref in guesses:
        try:
            v,l=advance(model,u,h,T,reference_guess=ref,vacuum_scaled=True,observer=observe)
            return v,l["residual"],1,len(errors)
        except (RuntimeError,FloatingPointError,RuntimeWarning) as e:
            errors.append(e)
    if continuation and depth==0:
        # These trial roots all solve backward Euler for the same unchanged
        # start moments. Only the final full-h solution is applied.
        theta=None
        try:
            for factor in (1/64,1/16,1/4,1.):
                theta,l=advance(model,u,h*factor,T,reference_guess=theta is None,
                                vacuum_scaled=True,initial_theta=theta)
            return theta,l["residual"],1,len(errors)
        except (RuntimeError,FloatingPointError,RuntimeWarning) as e:
            errors.append(e)
    if depth==5: raise errors[-1]
    v,r1,n1,f1=collide(model,u,h/2,T,depth+1,False,verbose)
    w,r2,n2,f2=collide(model,v,h/2,T,depth+1,False,verbose)
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

def run(eos,row,snap,nx=12,stop=410.,prompt_until=8.,orders=(3,2,2,8),off=False,radiation_table=None,dark_gravity=False,refinement=0,preroll_stop=None,preroll_steps=0,compact=False,use_pivot=False,prompt_only=False,progress=False,continuation=False,verbose=False,crop_edge=None):
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
    dropped=None
    if crop_edge is not None:
        target=crop_edge*row["initial_visible_T_GeV"]/row["final_a"]
        edge_index=int(np.argmin(abs(q-target)))
        assert edge_index>0 and abs(q[edge_index]-target)<2e-13
        old=Late(Collision(q),H(1.),[],False)
        counts_full,E_full=totals(old.k,old.moment(pack(th))[0],1.)
        assert max(abs(counts_full/np.asarray(base["species_counts_GeV3"])-1.))<4e-12
        q=q[:edge_index+1]
        th=[t[:,:edge_index].copy() for t in th]
        now=Late(Collision(q),H(1.),[],False)
        counts_keep,E_keep=totals(now.k,now.moment(pack(th))[0],1.)
        dropped=dict(q_over_TR=crop_edge,count_fraction=(1.-counts_keep/counts_full).tolist(),
                     a4_energy_fraction=(1.-E_keep/E_full).tolist())
    original_q=q.copy()
    q,th=subdivide(q,th,refinement)
    ks=Collision(q); model=Late(ks,H(1.),[],False)
    N0,E0=totals(ks,model.moment(pack(th))[0],1.)
    expected=np.asarray(base["species_counts_GeV3"])
    if dropped is not None: expected*=1.-np.asarray(dropped["count_fraction"])
    assert max(abs(N0/expected-1.))<4e-12
    init,_=report(model,th,1.,N0)
    # Use identical cold cutoff after subdivision; sum whole exact bins.
    def cold_fraction(mod,t):
        U=mod.moment(pack(t))[0].reshape(2,mod.n,2)
        return float(U[0,q[1:]<=initial_first_boundary*(1+1e-12),0].sum()/U[0,:,0].sum())
    init["scalar_count_below_archived_first_edge"]=cold_fraction(model,th)
    edge=16.*row["initial_visible_T_GeV"]/row["final_a"]
    idx=q[:-1]>=edge*(1-1e-12)
    def exterior(mod,t):
        N=mod.moment(pack(t))[0].reshape(2,mod.n,2)[:,:,0]
        return (N[:,idx].sum(axis=1)/N.sum(axis=1)).tolist()
    init["exterior_count_fractions_beyond_q16"]=exterior(model,th)
    e0=E0.sum()-np.array([mp,m])@N0
    reference_energy=(base["retained_dark_kinetic_energy_GeV4"] if dropped is None else
                      E_keep.sum()-np.array([mp,m])@counts_keep)
    assert abs(e0/reference_energy-1.)<4e-11
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
            if not prompt_only:
                setup(max(T,base["common_T_GeV"]/b),kk,T,first.eta,*orders,32.,.82,cache=panels,cache_only=compact)
            co=Late(kk,hb,panels,use_prompt,compact,use_pivot)
            uu,r,acc,rej=collide(co,u,y-x,T,continuation=continuation,verbose=verbose)
            maxres=max(maxres,r);accepted+=acc;rejected+=rej
            dE=totals(kk,co.moment(uu)[0],b)[1]-totals(kk,co.moment(u)[0],b)[1]
            Q4+=dE;Q3+=dE/b
            th=co.arrays(uu)
        th,t,*_=redshift(q,b,c,th);pressure+=t
        k=Collision(q/c)
        diagnostic=Late(k,H(c),[],False)
        summary,E=report(diagnostic,th,c,N0)
        summary["scalar_count_below_archived_first_edge"]=cold_fraction(diagnostic,th)
        summary["exterior_count_fractions_beyond_q16"]=exterior(diagnostic,th)
        summary.update(H_GeV=float(H(c,E.sum()/c**4)),
                       photon_temperature_GeV=tab.state(c)[0] if tab else base["actual_visible_T_GeV"]/c)
        steps.append(summary)
        if progress:
            import sys
            print(f"b={c:.9g} accepted={accepted} rejected={rejected}",file=sys.stderr,flush=True)
    ledger=max(abs(E-E0-pressure-Q4))/E0.sum()
    balance=abs(Q4.sum())/E0.sum()
    assert max(ledger,balance,max(x["count_defect"] for x in steps))<5e-9
    assert np.min(np.diff([init["dark_entropy_a3"]]+[x["dark_entropy_a3"] for x in steps]))>=-5e-11*init["dark_entropy_a3"]
    return dict(nsteps=nx,elastic_orders=orders,elastic=not off,prompt_until=prompt_until,
                initial=init,steps=steps,
                explicitly_discarded_archive_exterior=dropped,
                compact_elastic=compact,pivoted_prompt=use_pivot,prompt_only=prompt_only,residual_continuation=continuation,
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
    p.add_argument("--compact-elastic",action="store_true")
    p.add_argument("--pivot-prompt",action="store_true")
    p.add_argument("--prompt-only",action="store_true")
    p.add_argument("--progress",action="store_true")
    p.add_argument("--continuation",action="store_true")
    p.add_argument("--newton-log",action="store_true")
    p.add_argument("--crop-edge",type=float,help="diagnostic: retain only through this existing q/TR edge")
    a=p.parse_args()
    b=json.loads(Path(a.run).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    print(json.dumps(run(StandardModel(a.eos),b,s,a.steps,a.stop,a.prompt_until,a.orders,a.off,a.radiation_table,a.dark_gravity,a.refinement_stride,a.preroll_stop,a.preroll_steps,a.compact_elastic,a.pivot_prompt,a.prompt_only,a.progress,a.continuation,a.newton_log,a.crop_edge),indent=2,allow_nan=False))
