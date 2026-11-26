"""Late finite-shell check of the thermal-moment evolution, below pole only."""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import from_row,fugacity
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_combined_centered_implicit_reference import Combined
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import m,mp


class ElasticOnly(Combined):
    def reference(self,q,guess):
        U,_,_=self.moment(q)
        n=U.reshape(2,self.n,2)[:,:,0].sum(axis=1)
        self.t=guess
        # An arbitrary common-temperature detailed-balance reference
        # suffices for algebraic centering, even if the actual species
        # are far apart in temperature. Start Newton from the actual state.
        from continuum_dark_normal_capacity_reference import critical_T
        self.t=max(self.t,1.01*critical_T(n[0],20))
        et=[fugacity(self.t,C.mass,n[j],j==0,20) for j,C in enumerate(self.k.c)]
        self.eta=[C.mass/self.t+et[j] for j,C in enumerate(self.k.c)]
        ref=pack([np.array([C.low/self.t-et[j],C.dk/self.t]) for j,C in enumerate(self.k.c)])
        return ref,None,self.t

    def centered(self,z,ref,R0):
        return self.elastic(z)


def totals(kin,U,a):
    nc=len(kin.c[0].vol)
    v=U.reshape(2,nc,2)
    N=v[:,:,0].sum(axis=1)*a**3
    E=np.array([(v[j,:,0]*(C.low+C.mass)+v[j,:,1]*C.dk).sum() for j,C in enumerate(kin.c)])*a**4
    return N,E


def run(eos,base,two_temp,start_req=18.,stop=410.,steps=10,mesh=16,elastic=True,orders=(3,2,2,8),tail=32.):
    c=from_row(base,eos,20)
    ref=min(two_temp["steps"],key=lambda q:abs(np.log(q["b"]/start_req)))
    start,T0=ref["b"],ref["T_GeV"]
    assert start<stop
    N_target=np.asarray(c["species_counts_GeV3"])
    q=.05/base["final_a"]*np.r_[0.,np.geomspace(.005,16.,mesh)]
    k=Collision(q/start)
    H=lambda a:np.sqrt(np.pi**2/45)*(c["actual_visible_T_GeV"]/a)**2/MPL
    th=[]
    for i,C in enumerate(k.c):
        et=fugacity(T0[i],C.mass,N_target[i]/start**3,i==0,20)
        th.append(np.array([C.low/T0[i]-et,C.dk/T0[i]]))
    ent=ElasticOnly(k,H(start),[])
    U,_,_=ent.moment(pack(th))
    N0,E0=totals(k,U,start)
    pressure=np.zeros(2); energy_exchange=np.zeros(2); heat_a3=np.zeros(2)
    max_res=0.; nsum=0
    ncol=0
    xs=np.linspace(np.log(start),np.log(stop),steps+1)
    for x,y in zip(xs,xs[1:]):
        a,b,d=np.exp([x,(x+y)/2,y])
        th,pr,*_=redshift(q,a,b,th);pressure+=pr
        kk=Collision(q/b)
        if elastic:
            first=ElasticOnly(kk,H(b),[])
            X=pack(th)
            ref0,_,refT=first.reference(X,max(T0)*start/b)
            af=first.arrays(ref0)
            eta=[(C.mass+C.low[0])/refT-af[j][0,0] for j,C in enumerate(kk.c)]
            panels=[]
            setup(max(refT,max(T0)*start/b),kk,refT,eta,*orders,tail,.82,cache=panels)
            ncol+=len(panels)
            co=ElasticOnly(kk,H(b),panels)
            ux,_,_=co.moment(X)
            u1,stat=advance(co,X,y-x,refT,reference_guess=False,vacuum_scaled=True)
            ua,_,_=co.moment(u1)
            dE=totals(kk,ua,b)[1]-totals(kk,ux,b)[1]
            energy_exchange+=dE
            heat_a3+=dE/b
            th=co.arrays(u1)
            max_res=max(max_res,stat["residual"]); nsum+=stat["newton_iterations"]
        th,pr,*_=redshift(q,b,d,th);pressure+=pr
    kk=Collision(q/stop); co=ElasticOnly(kk,H(stop),[])
    ux,_,_=co.moment(pack(th))
    N1,E1=totals(kk,ux,stop)
    volumes=ux.reshape(2,-1,2)
    p=q/stop
    # First bin is [0,.005] in initial T_R units.
    low=volumes[0,0,0]/volumes[0,:,0].sum()
    central=abs(E1-E0-pressure-energy_exchange)/E0.sum()
    cn=max(abs(N1/N0-1))
    assert max(central.max(),cn,abs(energy_exchange.sum())/E0.sum())<7e-9
    from late_two_temperature_expansion_reference import inverse
    kf=(E1/stop-np.array([mp,m])*N1)/stop**3
    temps=[inverse(N1[j]/stop**3,kf[j],j,"K",T0[j]*start/stop)["T"] for j in range(2)]
    from continuum_dark_normal_capacity_reference import moments
    sratio=(N1[0]/stop**3)/moments(temps[0],mp,0.,True)[0]
    initK=E0.sum()/start-np.array([mp,m])@N0
    return dict(start=start,stop=stop,mesh=mesh,steps=steps,quadrature=orders,elastic=elastic,
                initial_moment_equivalent_T_GeV=T0,
                scalar_count_over_saturation_at_moment_equivalent_T=float(sratio),
                scalar_a3_heat_over_initial_total_a3_kinetic_energy=float(heat_a3[0]/initK),
                minimum_resolved_scalar_entropy_exponent=float(np.min(th[0][0]+np.minimum(0.,th[0][1]))),
                initial_shell_count_loss=(1-N0/N_target).tolist(),equivalent_final_T_GeV=temps,
                final_scalar_bin_fraction_below_q_over_TR_005=float(low),collective_pressure_energy_defect=central.tolist(),
                represented_count_defect=float(cn),signed_collision_energy_over_initial_total_E=(energy_exchange/E0.sum()).tolist(),
                max_solver_residual=max_res,newton_iterations=nsum,geometric_panels=ncol)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("base");p.add_argument("twotemp")
    p.add_argument("--steps",type=int,default=10);p.add_argument("--mesh",type=int,default=16);p.add_argument("--off",action="store_true")
    p.add_argument("--orders",type=int,nargs=4,default=(3,2,2,8))
    p.add_argument("--start",type=float,default=18.);p.add_argument("--stop",type=float,default=410.)
    a=p.parse_args()
    b=json.loads(Path(a.base).read_text())["rows"][-1]
    t=json.loads(Path(a.twotemp).read_text())["rows"][-1]
    print(json.dumps(run(StandardModel(a.eos),b,t,a.start,a.stop,a.steps,a.mesh,not a.off,a.orders),indent=2,allow_nan=False))
