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

def run(row,eos,nx,extent=.07,scale=2.3e-6,centered=False,warm=False,elastic=False):
    assert not elastic or centered
    vis,H,k,N,t,Q,th=configuration(row,eos,8,6.5)
    q=k.c[0].p
    omitted=[0.,0.]
    if warm:
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
    E0=e_species(k,u0,1.)
    N0=u0.reshape(2,cnt,2)[:,:,0].sum(axis=1)
    dE=np.zeros(2);Qadd=np.zeros(2)
    dSadv=dScol=0.;er=0.;en=0.;entropy_r=0.;nonlinear=0
    dx=extent/nx
    s0=eos.state(vis)[2]
    for i in range(nx):
        a,b,c=np.exp(dx*np.array([i,i+.5,i+1.]))
        th,g,x,y,s=redshift(q,a,b,th)
        dE+=g;er=max(er,x);en=max(en,y);dSadv+=s
        kin=Collision(q/b,order=5)
        v=brentq(lambda t:eos.state(t)[2]/s0-b**-3,vis*.7,vis)
        if centered:
            from joint_resonant_centered_stiff_reference import Centered
            m0=Centered(kin,H)
        else:m0=Model(kin,H)
        v0,s1,_=m0.moment(pack(th))
        Et=e_species(kin,v0,b)
        h=np.sqrt((eos.state(v)[0]+Et.sum()/b**4)/3)/MPL
        m0.H=h
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
                w,log=advance(m0,pack(th),dx*scale,v,reference_guess=warm)
            else:
                w,log=step(m0,pack(th),dx*scale)
            nonlinear+=log["newton_iterations"]
            th=m0.arrays(w)
        v1,s2,_=m0.moment(pack(th))
        Qadd+=e_species(kin,v1,b)-Et
        dScol+=(s2-s1)*b**3
        th,g,x,y,s=redshift(q,b,c,th);dE+=g;er=max(er,x);en=max(en,y);dSadv+=s
    af=np.exp(extent)
    kn=Collision(q/af,order=5)
    m1=Model(kn,H)
    U,S,_=m1.moment(pack(th))
    Ef=e_species(kn,U,af)
    n1=U.reshape(2,cnt,2)[:,:,0].sum(axis=1)*af**3
    errE=np.max(abs(Ef-E0-dE-Qadd))/E0.sum()
    errQ=abs(sum(Qadd))/E0.sum()
    errS=abs((S*af**3-S0)-dScol-dSadv)/S0
    assert max(errE,errQ,errS,er,np.max(abs(n1/N0-1)))<2e-9
    Tin,eq=maxent(kn.c,n1*af**-3,Ef.sum()*af**-4,vis/af)
    rem=[float(.5*np.sum(abs(U[2*j*cnt:2*(j+1)*cnt:2]-eq[j][2]))/(n1[j]*af**-3)) for j in (0,1)]
    return dict(nsteps=nx,final_a=af,common_equilibrium_T=Tin,centered_residual=centered,
      final_number_shell_total_variation_from_closed_box_equilibrium=rem,
      omitted_source_number_fractions=omitted,comoving_support_over_initial_visible_T=(q[[0,-1]]/vis).tolist(),
      initialization="warm source-bin counts" if warm else "manufactured near equilibrium",below_pole_elastic=elastic,manufactured_collision_amplitude_relative_to_actual=scale,
      total_number_defect=float(np.max(abs(n1/N0-1))),split_energy_ledger_defect=float(errE),
      collision_opposite_heat_defect=float(errQ),integrated_entropy_ledger_defect=float(errS),
      comoving_q_scalar_from_collisions_over_initial_energy=Qadd[0]/E0.sum(),
      redshift_trace_contribution_over_initial_energy=(dE/E0.sum()).tolist(),
      projection_entropy_change_relative_to_initial=dSadv/S0,
      collision_entropy_change_relative_to_initial=dScol/S0,
      maximum_moment_reconstruction_defect=en,nonlinear_iterations=nonlinear)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    p.add_argument("--physical-amplitude",action="store_true");p.add_argument("--warm-source-init",action="store_true");p.add_argument("--with-elastic",action="store_true")
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    rows=([run(r,eos,3,scale=1.,centered=True,warm=a.warm_source_init,elastic=a.with_elastic),run(r,eos,6,scale=1.,centered=True,warm=a.warm_source_init,elastic=a.with_elastic)]
          if a.physical_amplitude else [run(r,eos,3,scale=0.),run(r,eos,3),run(r,eos,6)])
    print(json.dumps(dict(rows=rows,
       scope=(("Warm source-bin initial projection with both ends cut. " if a.warm_source_init else "Manufactured near-equilibrium shape. ")+("Prompt return and below-pole elasticity; " if a.with_elastic else "Prompt return only; ")+"injection disabled.")),
       indent=2,allow_nan=False))
