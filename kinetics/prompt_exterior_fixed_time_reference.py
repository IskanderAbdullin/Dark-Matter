"""Short fixed-background transient of the resolved on-shell exterior."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from continuum_dark_normal_capacity_reference import from_row,solve,fugacity
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import TabulatedRadiation
from late_soft_two_moment_test_reference import totals
from pivoted_prompt_reference import Pivoted
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import m,mp


def run(eos,table,row,shape,extent,rtol=5e-9,method="Radau"):
    cap=from_row(row,eos,20)
    T=cap["actual_visible_T_GeV"]
    Ekin=cap["retained_dark_kinetic_energy_GeV4"]
    n=np.asarray(cap["species_counts_GeV3"])
    rh=np.array([mp,m])@n+Ekin
    H=np.sqrt((TabulatedRadiation(table,T).state(1.)[1]+rh)/3)/MPL
    q=np.asarray(shape["lab_momentum_edges_GeV"])
    boundary=16.*row["initial_visible_T_GeV"]/row["final_a"]
    mask=q[:-1]>=boundary*(1-1e-12)
    assert np.any(abs(q-boundary)<5e-13)
    model=Pivoted(Collision(q),H)
    th=[np.array(shape["entropy_variable_"+j]) for j in ("scalar","fermion")]
    X=pack(th)
    U,S0,_=model.moment(X)
    n0,E0=totals(model.k,U,1.)
    common=solve(n,Ekin)["common_T_GeV"]
    et=[fugacity(common,ma,n[j],j==0,20) for j,ma in enumerate((mp,m))]
    ref=pack([np.array([c.low/common-et[j],c.dk/common]) for j,c in enumerate(model.k.c)])
    R0=np.exp(sum(et)+((mp+m)-model.k.E[model.k.y])/common)
    Z=X-ref
    def fold(C,F):
        a,b,d=C.T
        den=a*d-b*b
        assert min(den)>0.
        return (-np.stack([d*F[:,0]-b*F[:,1],a*F[:,1]-b*F[:,0]],axis=1)/den[:,None]).ravel()
    def rhs(t,z):
        F=model.centered(z,ref,R0)[0].reshape(-1,2)
        _,_,C=model.moment(ref+z)
        return fold(C,F)
    def jac(t,z):
        F,jv=model.centered(z,ref,R0)
        _,_,C=model.moment(ref+z)
        u=fold(C,F.reshape(-1,2)).reshape(-1,2)
        h=[]
        for j,d in enumerate(model.k.c):
            theta=model.arrays(ref+z)[j]
            s=theta[0,:,None]+theta[1,:,None]*d.u
            if j==0:
                f=np.exp(-s)/(-np.expm1(-s)); sign=1
            else:
                from scipy.special import expit
                f=expit(-s); sign=-1
            w=d.w*d.vol[:,None]*f*(1+sign*f)*(1+sign*2*f)
            h.append(np.stack([np.sum(w*d.u**i,axis=1) for i in range(4)],axis=1))
        h0,h1,h2,h3=np.concatenate(h).T
        def direction(d):
            da,db=d.reshape(-1,2).T
            a,b,c=-(h0*da+h1*db),-(h1*da+h2*db),-(h2*da+h3*db)
            extra=np.stack([a*u[:,0]+b*u[:,1],b*u[:,0]+c*u[:,1]],axis=1)
            return fold(C,jv(d).reshape(-1,2)+extra)
        return np.stack([direction(x) for x in np.eye(len(z))],axis=1)
    times=extent*np.r_[0.,np.geomspace(1e-3,1.,16)]
    if method=="Radau":
        a=np.cos(0.37*np.arange(len(Z)))
        eps=3e-6
        fd=(rhs(0.,Z+eps*a)-rhs(0.,Z-eps*a))/(2*eps)
        jac_error=np.linalg.norm(fd-jac(0.,Z)@a)/np.linalg.norm(fd)
        assert jac_error<4e-7
        extra=dict(jac=jac)
    else:
        jac_error=None
        extra={}
    fit=solve_ivp(rhs,(0.,extent),Z,method=method,rtol=rtol,atol=rtol*1e-2,t_eval=times,**extra)
    assert fit.success
    rows=[]
    for t,z in zip(fit.t,fit.y.T):
        Q,S,_=model.moment(ref+z)
        N,E=totals(model.k,Q,1.)
        F=model.centered(z,ref,R0)[0]
        resh=Q.reshape(2,-1,2)
        flow=F.reshape(2,-1,2)
        tot=resh[:,mask,0].sum(axis=1)/N
        cold=resh[0,0,0]/N[0]
        rows.append(dict(hubble_time=float(t),exterior_fractions=tot.tolist(),
                    net_instantaneous_exterior_count_over_NH=(flow[:,mask,0].sum(axis=1)/N).tolist(),
                    scalar_first_cell_fraction=float(cold),
                    max_number_relative_defect=float(max(abs(N/n0-1.))),
                    total_energy_relative_defect=float(abs(E.sum()/E0.sum()-1.)),
                    relative_entropy_gain=float(S/S0-1.)))
    assert max(x["max_number_relative_defect"] for x in rows)<2e-8
    assert max(x["total_energy_relative_defect"] for x in rows)<2e-8
    assert min(np.diff([x["relative_entropy_gain"] for x in rows]))>-3e-11
    return dict(rows=rows,rtol=rtol,solver_evaluations=fit.nfev,method=method,analytic_rhs_JVP_error=jac_error)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("table");p.add_argument("base");p.add_argument("shape")
    p.add_argument("--method",choices=["Radau","DOP853"],default="Radau")
    p.add_argument("--extent",type=float,default=3e-8);p.add_argument("--rtol",type=float,default=5e-9)
    a=p.parse_args();row=json.loads(Path(a.base).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    print(json.dumps(run(StandardModel(a.eos),a.table,row,s,a.extent,a.rtol,a.method),indent=2,allow_nan=False))
