"""Bound the instantaneous warm archived stress energy, including exterior tails."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import snapshot_functions
from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import m,mp

MASSES=np.array([mp,m])


def retained(s,n,b=1.):
    e=np.asarray(s["lab_momentum_edges_GeV"])
    x,w=roots_legendre(n)
    p=(e[1:]+e[:-1])[:,None]/2+np.diff(e)[:,None]*x/2
    wt=np.diff(e)[:,None]*w/2*p*p/(2*np.pi**2)
    f,_,g,_=snapshot_functions(s)
    result=[]
    for mass,occ,spin in zip(MASSES,(f,g),(1,4)):
        h=np.hypot(p,mass)
        wt2=spin*wt*occ(h)/b**3
        pr=p/b
        er=np.hypot(pr,mass)
        r=float(np.sum(wt2*er))
        p0=float(np.sum(wt2*pr*pr/(3*er)))
        num=float(wt2.sum())
        result.append((num,r,p0))
    a=np.asarray(result)
    assert np.max(abs(a[:,0]*(s["scale_factor"]*b)**3/np.asarray(s["represented_species_count_a3"])-1.))<5e-11
    return a


def from_ledger(l,a4,number,a):
    rest=float(MASSES@number/a**3)
    lo=l["exterior_a4_energy_one_mean_per_cohort"]/a4
    hi=l["exterior_a4_energy_distribution_independent_upper_bound"]/a4
    assert 0<=rest<=lo<=hi
    return a,rest,lo-rest,hi-rest


def history(ledger,at):
    E=np.zeros(2)
    K=np.zeros(2)
    rest=0.
    for a,M,lo,hi in ledger:
        b=at/a
        assert b>=1
        rest+=M/b**3
        K+=np.array([lo/b**5,hi/b**4])
    E=rest+K
    # p^2/E_total = E_kin(1+m/E_total), bounded between
    # E_kin and 2 E_kin; furthermore P<=rho/3 for positive particles.
    P=np.array([K[0]/3,min(2*K[1]/3,E[1]/3)])
    assert np.all(E>=0) and 0<=P[0]<=P[1]
    return E,P


def test_envelope():
    rng=np.random.default_rng(1624)
    p=np.exp(rng.uniform(-22.,0.,(2,4000)))
    n=rng.uniform(size=p.shape)*1e-13
    base=np.hypot(p,MASSES[:,None])
    rest=float(np.sum(n*MASSES[:,None]))
    k=float(np.sum(n*p*p/(base+MASSES[:,None])))
    tight=[]
    for b in (1.,1.1,3.,10.,100.,1000.):
        q=p/b
        E=np.hypot(q,MASSES[:,None])
        exact=float(np.sum(n*E)/b**3)
        P=float(np.sum(n*q*q/(3*E))/b**3)
        x,y=history([(1.,rest,k,k)],b)
        assert x[0]<=exact*(1+3e-15) and exact<=x[1]*(1+3e-15)
        assert y[0]<=P*(1+3e-15) and P<=y[1]*(1+3e-15)
        tight.append(dict(expansion=b,energy_interval_over_exact=(x/exact).tolist(),
             pressure_interval_over_exact=(y/P).tolist()))
    return tight


def summary(s,ledger,EOS,b=1.,T=None):
    T=s["T_visible_GeV"] if T is None else T
    n=retained(s,40,b)
    q=retained(s,72,b)
    assert np.max(abs(n/q-1))<4e-10
    ae=s["scale_factor"]*b
    e,p=history(ledger,ae)
    v=TabulatedRadiation(EOS,T).state(1.)[1]
    x=(q[:,1].sum()+e)/v
    y=(q[:,2].sum()+p)/v
    delta=x/(np.sqrt(1+x)+1)
    return dict(T_visible_GeV=T,absolute_scale_factor=ae,
                retained_species_order=["scalar","stable"],
                retained_number_energy_pressure_columns=["number [GeV^3]","energy [GeV^4]","pressure [GeV^4]"],
                retained_number_energy_pressure=[list(row) for row in q],
                exterior_energy_ratio_to_visible_interval=(e/v).tolist(),
                exterior_pressure_ratio_to_visible_interval=(p/v).tolist(),
                instantaneous_dark_to_visible_energy_ratio_interval=x.tolist(),
                instantaneous_dark_w_envelope=[float(y[0]/x[1]),float(y[1]/x[0])],
                deltaH_over_dark_free_radiation_at_same_T_interval=delta.tolist(),
                retained_quadrature_relative_change=float(np.max(abs(n/q-1.))))


def read(p):
    return json.loads(Path(p).read_text())


def run(EOS,prod,budget,shape,stages,branches,free_T):
    B=read(budget)
    b=(B["cases"] if isinstance(B,dict) else B)[-1]
    assert b["enabled"] and b["signed_visible_portal"]
    s=read(shape)["snapshots"][-1]
    a=s["scale_factor"]
    n=np.asarray(s["represented_species_count_a3"])
    c=np.asarray(b["retained_fixed_Omega_h2"])/(MASSES*n)
    assert abs(c[0]/c[1]-1)<3e-12
    n_ex=np.asarray(b["exterior_fixed_Omega_h2"])/(MASSES*c)
    start=read(prod)["rows"][-1]
    assert abs(s["T_visible_GeV"]/start["final_visible_temperature_GeV"]-1)<1e-12
    ledger=[from_ledger(start["visible_bath_feedback"],a**4,n_ex,a)]
    out=[summary(s,ledger,EOS)]
    def advance(s,ledger,r,v):
        p=r["endpoint"]
        assert abs(r["initial_temperature_GeV"]/s["T_visible_GeV"]-1)<1e-12
        assert r["steps"]==v["steps"]
        ar=v["scale_factor"]
        ncon=np.asarray(s["represented_species_count_a3"])
        new=ncon*(np.asarray(p["inclusive_counts_relative_to_start"])
                     -np.asarray(p["retained_counts_relative_to_start"]))
        assert abs(v["T_visible_GeV"]/p["photon_temperature_GeV"]-1)<1e-12
        assert abs(ar/s["scale_factor"]/p["b"]-1)<1e-12
        assert np.all(new>=0)
        local_a4=(ar/s["scale_factor"])**4
        ext=ledger+[from_ledger(r["exterior"],local_a4,new,ar)]
        return v,ext,summary(v,ext,EOS)
    for rfile,sfile in stages:
        s,ledger,point=advance(s,ledger,read(rfile)["rows"][-1],read(sfile)["snapshots"][-1])
        out.append(point)
    independent=[]
    for rfile,sfile in branches:
        R=read(rfile)["rows"]
        V=read(sfile)["snapshots"]
        assert len(R)==len(V)>0
        independent.append([dict(steps=r["steps"],**advance(s,ledger,r,v)[2]) for r,v in zip(R,V)])
    future=[]
    if free_T:
        assert branches
        t,ext,_=advance(s,ledger,R[-1],V[-1])
        bg=TabulatedRadiation(EOS,t["T_visible_GeV"])
        for T in free_T:
            assert np.exp(bg.lo)<T<t["T_visible_GeV"]
            be=float(np.exp((bg.s0-bg.entropy_log(np.log(T)))/3))
            assert abs(bg.state(be)[0]/T-1.)<1e-11
            future.append(summary(t,ext,EOS,be,T))
    return out,independent,future


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","production","budget","starting_shape"):
        p.add_argument(k)
    p.add_argument("--stage",action="append",nargs=2,metavar=("ROWS","SHAPE"),default=[])
    p.add_argument("--branch",action="append",nargs=2,metavar=("ROWS","SHAPE"),default=[])
    p.add_argument("--free-future",nargs="*",type=float,default=[],metavar="T_GEV")
    a=p.parse_args()
    result,branches,future=run(a.eos,a.production,a.budget,a.starting_shape,a.stage,a.branch,a.free_future)
    print(json.dumps(dict(epochs=result,independent_branches_from_last_epoch=branches,
          extra_collisionless_future_from_last_branch=future,
          independent_two_mass_redshift_test=test_envelope(),scope="Instantaneous homogeneous expansion source from represented stable and scalar archives, with collisionless earlier exterior cohorts bounded by their endpoint number/energy. Not a self-consistent re-evolution with all tails, a radiation-equivalent Delta Neff, a BBN prediction, or a statement about observational expansion constraints."),
          indent=2,allow_nan=False))
