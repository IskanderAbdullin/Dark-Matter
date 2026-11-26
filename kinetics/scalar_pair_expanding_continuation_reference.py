"""Pair-only conservative expanding continuation of finite-shell archives."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import gross, m, mp, snapshot_functions
from conditional_relic_normalization_reference import conversion, proxy
from frw_exterior_cohort_reference import Exterior
from joint_prompt_two_moment_evolution_reference import Cells
from joint_resonant_frw_split_reference import redshift
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL

MASS=np.array([mp,m])


def moments(q,theta):
    result=[]
    for j,mass in enumerate(MASS):
        c=Cells(q,mass,(1,4)[j])
        k=c.low[:,None]+c.dk[:,None]*c.u
        f=c.f(np.asarray(theta[j]),mass+k,np.arange(len(c.low))[:,None])
        result.append(np.array([np.sum(c.w*f,axis=1),np.sum(c.w*f*k,axis=1)])*c.vol)
    return np.array(result)


def reconstruct(q,U):
    return [Cells(q,mi,gi).reconstruct(u)[0] for mi,gi,u in zip(MASS,(1,4),U)]


def e_total(U):
    return float(np.sum(U[:,1])+MASS@U[:,0].sum(axis=1))


def rate(p,theta,order,reverse=True):
    f,_,g,_=snapshot_functions(dict(lab_momentum_edges_GeV=p,
             entropy_variable_scalar=theta[0],entropy_variable_fermion=theta[1]))
    _,Fg,F=gross(p,m,mp,g,f,False,*order,project=True)
    if reverse:
        _,Rg,R=gross(p,mp,m,f,g,True,*order,project=True)
    else:
        R=np.zeros_like(F); Rg=np.zeros_like(Fg)
    net=np.array([F[1]-R[0],R[1]-F[0]])
    assert np.all(net[:,:,[0,-1]]>=0.)
    N=net[:,0].sum(axis=1)
    scale=2*(Fg[1]+Rg[1])
    nerr=abs(N.sum())/max(scale,1e-280)
    eerr=abs(e_total(net))/max(Fg[2]+Rg[2],1e-280)
    assert max(nerr,eerr)<2e-12
    return net,max(nerr,eerr),Rg[1]/Fg[1]


def cold(q,b,theta,cut):
    from scipy.special import roots_legendre
    p=q/b
    if cut<=q[0]:
        return 0.
    x,w=roots_legendre(30)
    hi=np.minimum(p[1:],cut/b)
    use=hi>p[:-1]
    hi=hi[use]
    lo=p[:-1][use]
    f,_,_,_=snapshot_functions(dict(lab_momentum_edges_GeV=p,
                 entropy_variable_scalar=theta[0],entropy_variable_fermion=theta[1]))
    pu=((hi+lo)[:,None]+(hi-lo)[:,None]*x)/2
    return float(b**3*np.sum((hi-lo)[:,None]*w*pu*pu*f(np.hypot(pu,mp)))/(4*np.pi**2))


def run(s,table,n,extent,order=(3,7,6,2),pad=0.,fixed=False,reverse=True,scale=1.):
    T=s["T_visible_GeV"]
    clock=TabulatedRadiation(table,T)
    q=np.asarray(s["lab_momentum_edges_GeV"])
    tr=s["T_initial_GeV"]
    ar=s["scale_factor"]
    th=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    if pad:
        target=pad*tr/ar
        assert target>q[-1]
        extra=int(np.ceil(np.log(target/q[-1])/.34))
        q=np.r_[q,np.geomspace(q[-1],target,extra+1)[1:]]
        th=[np.concatenate([a,np.vstack([np.full(extra,1000.),np.zeros(extra)])],axis=1) for a in th]
    state=moments(q,th)
    U0=state.copy()
    e0=e_total(U0)
    n0=U0[:,0].sum(axis=1)
    entropy=2*np.pi**2/45*np.exp(clock.s0)
    coeff=conversion(table)["Omega_h2_per_GeV_times_yield"]
    trace=0.;rate_error=0.;z0max=0.
    ex=Exterior();exn=np.zeros((2,2))
    dx=extent/n
    for i in range(n):
        a,b,c=np.exp(dx*(i+np.array([0.,.5,1.]))) if not fixed else (1.,1.,1.)
        if a!=b:
            th,work,_,_,_=redshift(q,a,b,th);trace+=work.sum()
        ex.advance(b)
        p=q/b
        U=moments(p,th)
        H=np.sqrt(clock.state(b)[1]/3)/MPL
        if scale:
            rr,err,_=rate(p,th,order,reverse=reverse)
            trial=reconstruct(p,U+(dx*scale/H/2)*rr[:,:,1:-1])
            rr,er,z=rate(p,trial,order,reverse=reverse)
            rate_error=max(rate_error,err,er)
            z0max=max(z0max,z)
            delta=(dx*scale/H)*rr
            th=reconstruct(p,U+delta[:,:,1:-1])
            for j,mi in enumerate(MASS):
                for side,i_tail in enumerate((0,-1)):
                    dn=b**3*delta[j,0,i_tail]
                    dk=b**4*delta[j,1,i_tail]
                    exn[j,side]+=dn
                    ex.add(b,mi,j,side,p[i_tail],dn,dk)
        if b!=c:
            th,work,_,_,_=redshift(q,b,c,th);trace+=work.sum()
        ex.advance(c)
    Uf=moments(q/c,th)
    ef=e_total(Uf)*c**4
    exr=ex.result()
    nf=Uf[:,0].sum(axis=1)*c**3
    nall=nf+exn.sum(axis=1)
    enerr=(ef+exr["exterior_a4_energy_one_mean_per_cohort"]-trace-ex.trace-e0)/e0
    nerr=nall.sum()/n0.sum()-1.
    assert max(abs(enerr),abs(nerr))<5e-10
    cut=.005*tr/ar
    return dict(steps=n,hold_background_fixed=fixed,scale_factor=c,
                duration_in_fixed_Hubble_units=extent if fixed else None,
                initial_and_final_visible_temperatures_GeV=[T,clock.state(c)[0]],
                quadrature_order=order,reverse_enabled=reverse,
                included_pair_operator=bool(scale), support_q_over_TR_at_original_origin=(q[[0,-1]]*ar/tr).tolist(),
                initial_fixed_count_Omega_h2=proxy(n0/entropy,coeff),
                retained_fixed_count_Omega_h2=proxy(nf/entropy,coeff),
                inclusive_fixed_count_Omega_h2=proxy(nall/entropy,coeff),
                inclusive_species_counts_over_initial=(nall/n0).tolist(),
                retained_species_counts_over_initial=(nf/n0).tolist(),
                exterior_counts_over_initial_species=(exn/n0[:,None]).tolist(),
                scalar_kinetic_energy_per_retained_particle_GeV=Uf[0,1].sum()/Uf[0,0].sum(),
                scalar_cold_fraction_q_over_original_TR_lt_005=(cold(q,c,th,cut)/nf[0] if cut>q[0] else None),
                total_count_ledger_defect=nerr,energy_ledger_defect=enerr,
                max_instantaneous_pair_ledger_defect=rate_error,max_inverse_over_forward=z0max,
                exterior=exr)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table"); p.add_argument("snapshot")
    p.add_argument("--stop",type=float,default=1.2)
    p.add_argument("--fixed",action="store_true")
    p.add_argument("--steps",type=int,nargs="+",default=[6,12])
    p.add_argument("--quadrature",type=int,nargs=4,default=[3,7,6,2])
    p.add_argument("--pad",type=float,default=0.)
    p.add_argument("--no-inverse",action="store_true")
    p.add_argument("--no-pair",action="store_true")
    a=p.parse_args()
    extent=a.stop if a.fixed else np.log(a.stop)
    assert extent>0. and min(a.steps)>0 and min(a.quadrature)>1
    s=json.loads(Path(a.snapshot).read_text())["snapshots"][-1]
    states=[run(s,a.table,n,extent,tuple(a.quadrature),a.pad,a.fixed,
                not a.no_inverse,0. if a.no_pair else 1.) for n in a.steps]
    print(json.dumps(dict(states=states,scope="Prescribed radiation after input archive; only represented pair reactions and conservative momentum redshift. Exterior moment cohorts tracked collisionlessly, do not feed reactions. No earlier production, other dark or visible collisions, or general cosmological relic claim."),indent=2,allow_nan=False))
