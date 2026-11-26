"""FRW conservative redshift and implicit prompt collisions on controlled finite supports."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,roots_legendre
from joint_prompt_two_moment_evolution_reference import Cells,Collision,maxent
from joint_resonant_implicit_step_reference import Model,step
from joint_resonant_kinetic_linear_modes_reference import configuration
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def redshift(q,a,b,theta):
    # q is comoving; both species use fixed comoving boundaries.
    pr=[None,None];led=[];count_res=[];shapeerr=[];mixS=0.
    from three_species_decay_onsager_reference import mp,m
    ms=((mp,1),(m,4))
    de=np.log(b/a)
    dx,dw=roots_legendre(24)
    time=de*(1+dx)/2
    for j,(mass,g) in enumerate(ms):
        c=Cells(q/a,mass,g);f=Cells(q/b,mass,g)
        js=np.asarray(theta[j])
        w=c.w*c.vol[:,None]
        v=js[0,:,None]+js[1,:,None]*c.u
        ff=expit(-v) if g==4 else np.exp(-v)/(-np.expm1(-v))
        k0=c.low[:,None]+c.dk[:,None]*c.u
        p0=np.sqrt(k0*(k0+2*mass))
        p1=p0*a/b
        k1=p1*p1/(np.hypot(p1,mass)+mass)
        N=(w*ff).sum(axis=1)*(a/b)**3
        K=(w*ff*k1).sum(axis=1)*(a/b)**3
        ini=np.sum(w*ff*(k0+mass))*a**4
        late=np.sum(K+N*mass)*b**4
        integ=np.sum((w*ff)[:,:,None] * np.exp(time)*(mass*mass/np.hypot(mass,p0[:,:,None]*np.exp(-time)))
                    *(dw*de/2))*a**4
        assert abs((late-ini)-integ)<2e-10*max(ini,late)
        zz,err,ent=f.reconstruct([N,K])
        pr[j]=zz;led.append(float(integ));shapeerr.append(err)
        count_res.append(abs(N.sum()*b**3/(np.sum(w*ff)*a**3)-1))
        s0=np.sum(w*(v*ff-(np.log(-np.expm1(-v)) if g==1 else -np.logaddexp(0.,-v))))*a**3
        mix=ent*b**3-s0
        assert mix>=-1e-10*s0
        mixS+=mix
    return pr,np.array(led),max(count_res),max(shapeerr),mixS

def pack(th):return np.concatenate([np.asarray(x).T.ravel() for x in th])

def run(row,eos,nx,extent=.07,scale=2.3e-6,centered=False,warm=False,elastic=False,source=None,empty_start=0.,support_low=.15,empty_cells=18,vacuum_scaled=False,allow_empty_shells=False,include_zero_cell=False,snapshot=None,empty_reheat=.05,visible_feedback=False,trace_exterior=False,empty_stop=.9,empty_support_high=16.,empty_outer_bins=0):
    assert nx>0 and (centered or not elastic)
    if visible_feedback: assert empty_start and source is not None
    if not empty_start:
        vis,H,k,_,_,_,th=configuration(row,eos,8,6.5)
        q=k.c[0].p
    omitted=[0.,0.]
    astart=1.
    pulse_full=None
    first_pulse_empty_shells=None
    if empty_start:
        assert source is not None and not warm and 0<empty_start<extent
        vis=empty_reheat
        s0=eos.state(vis)[2]
        astart=np.exp(empty_start)
        tmid=brentq(lambda t:eos.state(t)[2]/s0-np.exp(-1.5*empty_start),min(.9,empty_stop)*vis,vis)
        H=np.sqrt(eos.state(tmid)[0]/3)/MPL
        if empty_outer_bins:
            assert empty_support_high>16.
            q=vis*np.r_[np.geomspace(support_low,16.,empty_cells+1),
                        np.geomspace(16.,empty_support_high,empty_outer_bins+1)[1:]]
        else:
            q=vis*np.geomspace(support_low,empty_support_high,empty_cells+1)
        if include_zero_cell:
            q=np.r_[0.,q]
        k=Collision(q/astart)
        # End-deposit one tiny physical first pulse. No positive occupation
        # floor: a numerically unresolved retained cell is rejected.
        pulse_full=(empty_start/H)*source.rates(tmid,q,astart)
        kept=np.stack([pulse_full[2,:,1:-1],pulse_full[:2,:,1:-1].sum(axis=0)])
        first_pulse_empty_shells=(kept[:,0,:]<=0.).sum(axis=1)
        if np.any(first_pulse_empty_shells):
            if not allow_empty_shells or min(kept[:,0,:].sum(axis=1))<=0.:
                raise ValueError("initial physical source misses a retained shell")
            # Zero is passed to Cells.reconstruct unchanged. Its finite
            # exponent sentinel evaluates to exactly zero in float64.
            # The collision solve must use a positive *trial* (not input)
            # to permit subsequent gain to populate these bins.
        th=[C.reconstruct(x)[0] for C,x in zip(k.c,kept)]
    elif warm:
        p=np.array(row["scalar_differential_bins"])*vis
        lo=np.searchsorted(p,.15*vis,"left")
        hi=np.searchsorted(p,15.9*vis,"right")-1
        ix=np.unique(np.r_[np.arange(lo,hi+1,6),hi])
        q=p[ix]
        k=Collision(q)
        ss=eos.state(vis)[2]
        th=[]
        for j,name in enumerate(("scalar_normalized_counts","stable_normalized_counts")):
            n=np.array(row[name])
            kept=np.add.reduceat(n[:hi],ix[:-1])
            omitted[j]=1-float(kept.sum())
            y=row["stable_yields"][1] if j==0 else sum(row["stable_yields"])
            occ=kept*(ss*y)/k.c[j].vol
            assert min(occ)>0 and (j==0 or max(occ)<1)
            a=np.log1p(occ)-np.log(occ) if j==0 else np.log1p(-occ)-np.log(occ)
            th.append(np.array([a,np.zeros_like(a)]))
    else:
        for j in (0,1):
            th[j][0]+=0.11*np.sin(np.arange(len(th[j][0]))*.81+.7*j)
    init=Model(k,H)
    cnt=len(q)-1
    u0,S0,_=init.moment(pack(th))
    def e_species(k,u,a):
        return a**4*np.array([u[2*j*cnt:2*(j+1)*cnt:2]@(c.mass+c.low)+u[2*j*cnt+1:2*(j+1)*cnt:2]@c.dk
                             for j,c in enumerate(k.c)])
    E0=e_species(k,u0,astart)
    N0=u0.reshape(2,cnt,2)[:,:,0].sum(axis=1)*astart**3
    dE=np.zeros(2);Qadd=np.zeros(2)
    dSadv=dScol=0.;er=0.;en=0.;nonlinear=0;nonlinear_max=0.;nonlinear_first=None
    source_energy=np.zeros(2);source_n=np.zeros(2)
    tails=np.zeros((3,2,2)); source_total=np.zeros((3,2))
    species_masses=np.array([k.c[1].mass,k.c[1].mass,k.c[0].mass])[:,None]
    def full_channels(values,a):
        tot=values.copy()
        tot[:,1]+=tot[:,0]*species_masses
        return tot*np.array([a**3,a**4])[None,:,None]
    if empty_start:
        first=full_channels(pulse_full,astart)
        first_tails=first[:,:,[0,-1]]
        first_all=first.sum(axis=2)
        kept_first=(first_all-first_tails.sum(axis=2))
        check=np.stack([kept_first[2],kept_first[:2].sum(axis=0)])
        assert np.max(abs(np.stack([N0,E0],axis=1)/check-1))<4e-10
    exterior_trace=None
    if trace_exterior:
        from frw_exterior_cohort_reference import Exterior
        exterior_trace=Exterior()
    def apply_exterior(rate,a):
        if exterior_trace is None: return
        for j,mass in enumerate(species_masses[:,0]):
            for side,idx in enumerate([0,-1]):
                n,k=rate[j,:,idx]
                if min(n,k)<0.:
                    raise ValueError("one-cohort exterior control does not absorb negative net tail rates")
                exterior_trace.add(a,mass,j,side,q[idx]/a,a**3*n,a**4*k)
    bath=None
    if visible_feedback:
        from visible_frw_energy_ledger_reference import VisibleBath
        bath=VisibleBath(eos,vis)
        bath.advance(astart)
        bath.exchange(float(first_all[:,1].sum()))
        if exterior_trace is not None:
            exterior_trace.advance(astart)
            apply_exterior(pulse_full,astart)
    source_entropy=0.
    signed_source=bool(source is not None and hasattr(source,"rates_at_state"))
    source_abs=np.zeros((3,2,2))
    dx=(extent-empty_start)/nx
    s0=eos.state(vis)[2]
    Tmin=max(np.nextafter(np.exp(eos.lo),np.inf),vis*min(.7,.9*empty_stop))
    for i in range(nx):
        a,b,c=np.exp(empty_start+dx*np.array([i,i+.5,i+1.]))
        th,g,x,y,s=redshift(q,a,b,th)
        dE+=g;er=max(er,x);en=max(en,y);dSadv+=s
        kin=Collision(q/b,order=5)
        if bath is not None:
            bath.advance(b)
            v=bath.T
        else:
            v=brentq(lambda t:eos.state(t)[2]/s0-b**-3,Tmin,vis)
        if centered:
            from joint_resonant_centered_stiff_reference import Centered
            m0=Centered(kin,H)
        else:m0=Model(kin,H)
        v0,s1,_=m0.moment(pack(th))
        Et=e_species(kin,v0,b)
        if exterior_trace is not None: exterior_trace.advance(b)
        outside=exterior_trace.energy(b) if exterior_trace is not None else 0.
        h=np.sqrt((eos.state(v)[0]+(Et.sum()+outside)/b**4)/3)/MPL
        m0.H=h
        if source is not None:
            # For the baseline the two midpoint increments share the
            # prescribed vacuum rate.  The signed source is re-evaluated
            # before each kick at its current represented occupations.
            fixed=None if signed_source or bath is not None else source.rates(v,q,b)
            def kick(before):
                temp=bath.T if bath is not None else v
                signed=source.rates_at_state(temp,q,b,before)*(dx/(2*h)) if signed_source else None
                if signed is not None:
                    for j in range(2):
                        source_abs[:,j]+=full_channels(signed[:,j],b).sum(axis=2)
                    sr=signed[:,0]-signed[:,1]
                else: sr=(source.rates(temp,q,b) if fixed is None else fixed)*(dx/(2*h))
                sr_tot=full_channels(sr,b)
                tails[:]+=sr_tot[:,:,[0,-1]]
                source_total[:]+=sr_tot.sum(axis=2)
                if bath is not None:
                    bath.exchange(float(sr_tot[:,1].sum()))
                apply_exterior(sr,b)
                kept=np.stack([sr[2,:,1:-1],sr[:2,:,1:-1].sum(axis=0)])
                source_n[:]+=b**3*kept[:,0,:].sum(axis=1)
                source_energy[:]+=b**4*np.array([(x[1]+C.mass*x[0]).sum() for x,C in zip(kept,kin.c)])
                U,S,_=m0.moment(pack(before))
                u=U.reshape(2,cnt,2)
                states=[]
                for j,C in enumerate(kin.c):
                    n=u[j,:,0]
                    K=n*C.low+u[j,:,1]*C.dk
                    reconstructed,_,_=C.reconstruct([n+kept[j,0],K+kept[j,1]])
                    states.append(reconstructed)
                _,after,_=m0.moment(pack(states))
                return states,(after-S)*b**3
            th,ds=kick(th);source_entropy+=ds
            if bath is not None:
                v=bath.T
            v0,s1,_=m0.moment(pack(th))
            Et=e_species(kin,v0,b)
        if elastic:
            from joint_below_pole_elastic_linear_modes_reference import setup
            from joint_combined_centered_implicit_reference import Combined
            ref,_,temp=m0.reference(pack(th),v)
            a0=m0.arrays(ref)
            eta=[(c.mass+c.low[0])/temp-a0[j][0,0] for j,c in enumerate(kin.c)]
            panels=[]
            setup(v,kin,temp,eta,4,3,2,10,30.,.82,cache=panels)
            m0=Combined(kin,h,panels)
        if scale:
            if centered:
                from joint_resonant_centered_stiff_reference import advance
                w,log=advance(m0,pack(th),dx*scale,v,
                              reference_guess=vacuum_scaled or warm or bool(empty_start and i>0 and support_low>=.04),
                              vacuum_scaled=vacuum_scaled)
                nonlinear_max=max(nonlinear_max,log["residual"])
                if i==0: nonlinear_first=log
            else:
                w,log=step(m0,pack(th),dx*scale)
            nonlinear+=log["newton_iterations"]
            th=m0.arrays(w)
        v1,s2,_=m0.moment(pack(th))
        Qadd+=e_species(kin,v1,b)-Et
        dScol+=(s2-s1)*b**3
        if source is not None:
            th,ds=kick(th);source_entropy+=ds
        th,g,x,y,s=redshift(q,b,c,th);dE+=g;er=max(er,x);en=max(en,y);dSadv+=s
        if bath is not None: bath.advance(c)
        if exterior_trace is not None: exterior_trace.advance(c)
    af=np.exp(extent)
    kn=Collision(q/af,order=5)
    m1=Model(kn,H)
    U,S,_=m1.moment(pack(th))
    Ef=e_species(kn,U,af)
    tf_base=brentq(lambda t:eos.state(t)[2]/s0-af**-3,Tmin,vis)
    tf=tf_base if bath is None else bath.T
    if snapshot is not None:
        snapshot.update(T_initial_GeV=vis, T_visible_GeV=tf, scale_factor=af,
                        lab_momentum_edges_GeV=(q/af).tolist(),
                        entropy_variable_scalar=np.asarray(th[0]).tolist(),
                        entropy_variable_fermion=np.asarray(th[1]).tolist(),
                        represented_species_count_a3=(U.reshape(2,cnt,2)[:,:,0].sum(axis=1)*af**3).tolist())
    # Independently reconcile in-domain inputs with the full rate less tails.
    if signed_source:
        assert np.max(abs(source_abs[:,0]-source_abs[:,1]-source_total)
                      /np.maximum(source_abs.sum(axis=1),1e-280))<3e-13
    kept_rates=source_total-tails.sum(axis=2)
    sums=np.stack([source_n,source_energy],axis=1)
    assert np.max(abs(sums-np.stack([kept_rates[2],kept_rates[:2].sum(axis=0)]))
                  /np.maximum(abs(sums),1e-100))<3e-13
    n1=U.reshape(2,cnt,2)[:,:,0].sum(axis=1)*af**3
    unit=E0.sum()+(source_energy.sum() if empty_start else 0.)
    errE=np.max(abs(Ef-E0-dE-Qadd-source_energy))/unit
    errQ=abs(sum(Qadd))/unit
    S0*=astart**3
    errS=abs((S*af**3-S0)-dScol-dSadv-source_entropy)/(S*af**3 if empty_start else S0)
    assert max(errE,errQ,errS,er,np.max(abs(n1/(N0+source_n)-1)))<2e-9
    visible_out=None
    if bath is not None:
        visible_out=bath.summary()
        exterior=first_tails[:,1].sum()+tails[:,1].sum()
        portal=first_all[:,1].sum()+source_total[:,1].sum()
        cross=abs(portal-bath.portal)/abs(portal)
        joint=((bath.energy-bath.E0-bath.trace)+Ef.sum()-dE.sum()+exterior)/bath.E0
        assert abs(joint)<4e-11 and cross<3e-13
        if exterior_trace is not None:
            ex=exterior_trace.result()
            assert abs(ex["exterior_a4_portal_injection"]-exterior)<3e-13*max(abs(exterior),1e-280)
            full=((bath.energy-bath.E0-bath.trace)+Ef.sum()-dE.sum()
                +ex["exterior_a4_energy_one_mean_per_cohort"]-ex["exterior_a4_pressure_trace"])/bath.E0
            assert abs(full)<4e-11
            ex["visible_dark_exterior_joint_ledger_defect_in_initial_visible_units"]=full
            visible_out.update(ex)
        visible_out.update(visible_baseline_T_same_a_GeV=tf_base,
                           relative_visible_temperature_shift=tf/tf_base-1.,
                           total_excluded_daughter_a4_energy=float(exterior),
                           total_portal_withdrawal_over_initial_visible_energy=portal/bath.E0,
                           combined_ledger_defect_in_initial_visible_units=joint,
                           visible_dark_portal_cross_defect_relative=cross)
    Tin,eq=maxent(kn.c,n1*af**-3,Ef.sum()*af**-4,vis/af)
    rem=[float(.5*np.sum(abs(U[2*j*cnt:2*(j+1)*cnt:2]-eq[j][2]))/(n1[j]*af**-3)) for j in (0,1)]
    soft_frac={}
    gx,gw=roots_legendre(32)
    ars=th[0]
    C=kn.c[0]
    for cut in (.005,.01,.025,.05,.1):
        if cut*vis<=q[0] or cut*vis>=q[-1]:
            continue
        edges=np.minimum(C.p[1:],cut*vis/af)
        valid=edges>C.p[:-1]
        lo=C.p[:-1][valid]
        hi=edges[valid]
        pp=(lo[:,None]*(1-gx)+hi[:,None]*(1+gx))/2
        kk=pp*pp/(np.hypot(pp,C.mass)+C.mass)
        zz=ars[0][valid,None]+ars[1][valid,None]*((kk-C.low[valid,None])/C.dk[valid,None])
        count=C.g*np.sum(gw*pp*pp*np.exp(-zz)/(-np.expm1(-zz))*(hi-lo)[:,None]/(4*np.pi**2))
        soft_frac[str(cut)]=float(count/(n1[0]/af**3))
    return dict(nsteps=nx,final_a=af,initial_visible_T_GeV=vis,common_equilibrium_T=Tin,centered_residual=centered,
      final_number_shell_total_variation_from_closed_box_equilibrium=rem,
      omitted_source_number_fractions=omitted,comoving_support_over_initial_visible_T=(q[[0,-1]]/vis).tolist(),
      initialization=("physical first source pulse; empty at a=1" if empty_start else ("warm source-bin counts" if warm else "manufactured near equilibrium")),below_pole_elastic=elastic,manufactured_collision_amplitude_relative_to_actual=scale,
      total_number_defect=float(np.max(abs(n1/(N0+source_n)-1))),split_energy_ledger_defect=float(errE),
      collision_opposite_heat_defect=float(errQ),integrated_entropy_ledger_defect=float(errS),
      comoving_q_scalar_from_collisions_over_initial_energy=Qadd[0]/E0.sum() if not empty_start else None,
      redshift_trace_contribution_over_initial_energy=(dE/E0.sum()).tolist() if not empty_start else None,
      projection_entropy_change_relative_to_initial=dSadv/S0 if not empty_start else None,
      collision_entropy_change_relative_to_initial=dScol/S0 if not empty_start else None,
      maximum_moment_reconstruction_defect=en,nonlinear_iterations=nonlinear,
      vacuum_scaled_residual=vacuum_scaled,
      maximum_accepted_centered_step_residual=nonlinear_max if centered and scale else None,
      first_collision_solver=nonlinear_first,
      represented_momentum_cells=cnt,
      innermost_positive_q_over_initial_T=q[1]/vis if include_zero_cell else q[0]/vis,
      final_resolved_scalar_count_fraction_below_q_over_initial_T=soft_frac,
      first_pulse_delta_log_a=empty_start,
      first_pulse_empty_shells=first_pulse_empty_shells.tolist() if empty_start else None,
      first_pulse_count_a3=N0.tolist() if empty_start else None,
      first_pulse_energy_a4=E0.tolist() if empty_start else None,
      first_pulse_count_and_total_energy_by_channel=first_all.tolist() if empty_start else None,
      first_pulse_count_and_total_energy_tails=first_tails.tolist() if empty_start else None,
      collision_q_scalar_over_scale=Qadd[0]/unit,
      scale="all kept source energy at injection" if empty_start else "initial total dark energy",
      final_scalar_over_dark_energy=Ef[0]/Ef.sum(),
      pressure_work_over_energy_scale=(dE/unit).tolist(),
      final_energy_over_scale=(Ef/unit).tolist(),
      final_counts_div_initial_visible_entropy=(n1/s0).tolist(),
      visible_injection=source is not None,
      final_visible_temperature_GeV=tf,
      visible_bath_feedback=visible_out,
      state_dependent_quantum_visible_source=signed_source,
      source_positive_negative_a3N_a4E=source_abs.tolist() if signed_source else None,
      comoving_energy_from_source_over_initial_total=(source_energy/E0.sum()).tolist() if not empty_start else None,
      comoving_count_from_source_over_initial_species=(source_n/N0).tolist() if not empty_start else None,
      comoving_source_channel_count_and_total_energy=source_total.tolist(),
      comoving_source_channel_count_and_total_energy_tails=tails.tolist())

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra',nargs='?')
    p.add_argument("--physical-amplitude",action="store_true");p.add_argument("--warm-source-init",action="store_true");p.add_argument("--with-elastic",action="store_true")
    p.add_argument("--visible-source",action="store_true");p.add_argument("--source-only",action="store_true")
    p.add_argument("--steps",type=int,nargs="+",default=[3,6])
    p.add_argument("--empty-pulse",type=float,default=0.)
    p.add_argument("--support-low",type=float,default=.15)
    p.add_argument("--empty-cells",type=int,default=18)
    p.add_argument("--empty-reheat",type=float,default=.05)
    p.add_argument("--empty-support-high",type=float,default=16.)
    p.add_argument("--empty-outer-bins",type=int,default=0)
    p.add_argument("--empty-stop",type=float,default=.9,help="reference visible final/initial temperature used to set final scale factor")
    p.add_argument("--vacuum-scaled-residual",action="store_true")
    p.add_argument("--allow-empty-shells",action="store_true")
    p.add_argument("--include-zero-cell",action="store_true")
    p.add_argument("--save-snapshots",help="write final finite-shell entropy variables for each run to this JSON file")
    p.add_argument("--quantum-visible-source",action="store_true")
    p.add_argument("--matched-vacuum-order",action="store_true")
    p.add_argument("--visible-backreaction",action="store_true")
    p.add_argument("--accurate-visible-eos",action="store_true")
    p.add_argument("--trace-exterior",action="store_true")
    a=p.parse_args();eos=StandardModel(a.eos)
    if a.empty_pulse: assert 0.<a.empty_stop<1.
    if a.trace_exterior: assert a.visible_backreaction
    if a.visible_backreaction or a.accurate_visible_eos:
        from visible_frw_energy_ledger_reference import ConservativeEOS
        eos=ConservativeEOS(eos)
    source=None
    if a.visible_source:
        from joint_visible_fd_source_reference import Visible
        source=Visible()
    assert a.empty_pulse>=0 and 0<a.support_low<a.empty_support_high and a.empty_cells>=4
    if a.include_zero_cell: assert a.allow_empty_shells
    if a.allow_empty_shells: assert a.vacuum_scaled_residual
    if a.vacuum_scaled_residual:
        assert a.physical_amplitude and a.empty_pulse and not a.source_only
    if a.matched_vacuum_order: assert a.visible_source and not a.quantum_visible_source
    if a.quantum_visible_source:
        assert a.visible_source and a.empty_pulse
        from joint_visible_quantum_source_reference import QuantumVisible
        source=QuantumVisible()
    elif a.matched_vacuum_order:
        from joint_visible_fd_source_reference import Visible
        source=Visible(order=(6,12,14),polar=16)
    if a.source_only: assert a.visible_source and not a.with_elastic
    if a.visible_source: assert a.physical_amplitude
    r=None
    if a.empty_pulse:
        assert a.visible_source and not a.warm_source_init
    else:
        if not a.spectra: p.error('spectra required except with --empty-pulse')
        r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    extent=np.log(eos.state(a.empty_reheat)[2]/eos.state(a.empty_stop*a.empty_reheat)[2])/3 if a.empty_pulse else .07
    final_shapes=[{} for _ in a.steps] if a.save_snapshots else [None for _ in a.steps]
    if a.save_snapshots: assert a.physical_amplitude and a.empty_pulse
    rows=([run(r,eos,n,extent=extent,empty_start=a.empty_pulse,support_low=a.support_low,empty_cells=a.empty_cells,vacuum_scaled=a.vacuum_scaled_residual,allow_empty_shells=a.allow_empty_shells,include_zero_cell=a.include_zero_cell,scale=0. if a.source_only else 1.,centered=True,warm=a.warm_source_init,elastic=a.with_elastic,source=source,snapshot=snapshot,empty_reheat=a.empty_reheat,visible_feedback=a.visible_backreaction,trace_exterior=a.trace_exterior,empty_stop=a.empty_stop,empty_support_high=a.empty_support_high,empty_outer_bins=a.empty_outer_bins) for n,snapshot in zip(a.steps,final_shapes)]
          if a.physical_amplitude else [run(r,eos,3,scale=0.),run(r,eos,3),run(r,eos,6)])
    if a.save_snapshots:
        for n,x in zip(a.steps,final_shapes): x.update(steps=n, below_pole_elastic=a.with_elastic)
        Path(a.save_snapshots).write_text(json.dumps(dict(snapshots=final_shapes),indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(rows=rows,
       scope=(("Warm source-bin initial projection with both ends cut. " if a.warm_source_init else ("Actually empty initial state on a finite fixed support. " if a.empty_pulse else "Manufactured near-equilibrium shape. "))+("Prompt return and below-pole elasticity; " if a.with_elastic else "Prompt return only; ")+(("state-dependent visible FD source enabled." if a.quantum_visible_source else "visible source enabled (vacuum final state).") if a.visible_source else "injection disabled."))),
       indent=2,allow_nan=False))
