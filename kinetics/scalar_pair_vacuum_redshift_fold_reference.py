"""Independent shape-frozen vacuum fold and optimistic all-history bound."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import gross, snapshot_functions, M, m, mp, y
from conditional_relic_normalization_reference import conversion, proxy
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import C, MPL


def fold(snapshot,table,stop,order,time_order):
    """Forward-only vacuum rate with collisionlessly redshifted shape."""
    T=snapshot["T_visible_GeV"]
    edges=np.array(snapshot["lab_momentum_edges_GeV"])
    _,_,f0,_=snapshot_functions(snapshot)
    clock=TabulatedRadiation(table,T)
    entropy=2*np.pi**2/45*np.exp(clock.s0)
    yields=np.array(snapshot["represented_species_count_a3"])/(snapshot["scale_factor"]**3*entropy)
    fee=conversion(table)["Omega_h2_per_GeV_times_yield"]
    bounds=np.unique(np.clip([1.,1.5,2.,3.,5.,10.,30.,60.,stop],1,stop))
    optical=[]
    x,w=roots_legendre(time_order)
    for left,right in zip(np.log(bounds[:-1]),np.log(bounds[1:])):
        integ=0.
        for u,g in zip((left+right)/2+(right-left)/2*x,(right-left)/2*w):
            b=np.exp(u)
            def f(E):
                p=np.sqrt(np.maximum((E-m)*(E+m),0.))
                return f0(np.hypot(m,b*p))
            n,ev,_=gross(edges/b,m,mp,f,np.zeros_like,False,*order,split_final=False)
            H=np.sqrt(clock.state(b)[1]/3)/MPL
            integ+=g*2*ev[1]/(n*H)
        optical.append(integ)
    tau=sum(optical)
    delta=yields[1]*tau/(1+tau)
    out=yields+np.array([delta,-delta])
    return dict(final_scale_factor=stop,order=order,time_order=time_order,
                scale_factor_panels=bounds.tolist(),optical_depths=optical,optical_depth=tau,
                counts_after_over_before=(out/yields).tolist(),
                initial_fixed_count_Omega_h2=proxy(yields,fee),
                final_fixed_count_Omega_h2=proxy(out,fee))


def optimistic(table,tr,full_fixed_stable_omega):
    """Allow the whole stable injection at reheating and maximal vacuum rate."""
    t,gs,_=np.loadtxt(table).T
    gmax=max(gs[:np.searchsorted(t,tr)+1])
    fee=conversion(table)["Omega_h2_per_GeV_times_yield"]
    Y=full_fixed_stable_omega/(fee*m)
    A=np.sqrt(2*C/3)/MPL
    bound=y**4/(8*np.pi*M*M)
    tau=bound*Y/2*(4*C/3)*gmax*tr/A
    return dict(T_R_GeV=tr,assumed_cumulative_stable_vacuum_Omega_h2=full_fixed_stable_omega,
                max_vacuum_optical_depth=tau,
                upper_increment_in_scalar_fixed_count_Omega_h2=full_fixed_stable_omega*mp/m*tau/(1+tau))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("snapshot")
    p.add_argument("--stops",type=float,nargs="+",default=[10.,60.])
    p.add_argument("--optimistic-total-stable",type=float)
    a=p.parse_args()
    assert min(a.stops)>1
    s=json.loads(Path(a.snapshot).read_text())["snapshots"][-1]
    rows=[]
    for b in a.stops:
        q=[fold(s,a.table,b,order,tm) for order,tm in (((3,8,8,14),5),((5,12,12,20),8))]
        assert abs(q[0]["optical_depth"]/q[1]["optical_depth"]-1.)<3e-5
        rows.extend(q)
    result=dict(rows=rows,scope="Prescribed radiation, no further source/collisions/inverse or stimulation; fixed redshifted fermion shape with uniformly evolving count. Not a relic prediction.")
    if a.optimistic_total_stable is not None:
        result["optimistic_vacuum_bound"]=optimistic(a.table,s["T_initial_GeV"],a.optimistic_total_stable)
    print(json.dumps(result,indent=2,allow_nan=False))
