"""Continuation of the dilute signed closure with expansion and direct visible input."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp
from warm_prompt_reheating_signed_slow_closure_reference import closed,signed


def visible(bg,T):
    u=np.log(T)
    rho=np.pi**2/30*T**4*np.exp(bg.gr(u))
    sent=2*np.pi**2/45*T**3*np.exp(bg.gs(u))
    return rho,sent


def source_table(Ts,n):
    a=[]
    ls=Source(((9,28,24),(14,46,32))[n],tail=72.,dark_mass=m)
    hs=Source(((9,28,24),(14,46,32))[n],tail=72.,dark_mass=M)
    for T in Ts:
        q,_,_,C,_=ls.event_moments(T)
        p,_,_,D,_=hs.event_moments(T)
        a.append([2*(C.sum()+D.sum()),2*D.sum(),(C*q).sum()+(D*p).sum()])
    return PchipInterpolator(np.log(Ts),np.log(a),axis=0)


def one(bg,row,q,closure_limit,include_source):
    T0=row["TF"];f0=np.log(bg.scale(T0))
    x_end=np.log(bg.scale(.035))-f0
    Y=np.array((sum(row["stable_yields"]),row["stable_yields"][1]))
    E=np.array(row["redshifted_means"])
    init=np.r_[Y,E[:2]@row["stable_yields"]+E[2]*Y[1]]
    source=source_table(np.geomspace(.0345,.0605,(21,39)[q]),q)
    order=(24,52)[q]
    def calc(x,u):
        T=np.exp(bg.v_to_u(f0+x))
        rv,sv=visible(bg,T)
        Td,N,z=closed(u[:2]*sv,u[2]*sv,T)
        H=np.sqrt((rv+u[2]*sv)/3)/MPL
        J=np.exp(source(np.log(T))) if include_source else np.zeros(3)
        W,_,_=signed(Td,z,order)
        rhs=np.r_[(W+J[:2])/(H*sv),J[2]/(H*sv)-3*Td*N.sum()/sv]
        pk=z*np.exp(-np.array((m,mp,M))/Td)
        return T,Td,N/sv,pk,rhs,rv,sv,np.array([J[2]/(H*sv),3*Td*N.sum()/sv])
    def boundary(x,u):
        return np.max(calc(x,u)[3])-closure_limit
    boundary.terminal=True
    boundary.direction=1
    sol=solve_ivp(lambda x,u:calc(x,u)[4],(0.,x_end),init,rtol=(1e-7,1e-9)[q],
                  atol=(1e-14,1e-16)[q],method="DOP853",max_step=.004,dense_output=True,events=boundary)
    assert sol.success and sol.y.min()>0
    # Verify the energy flow against an independent full quadrature on
    # the reconstructed trajectory, without subtracting two close totals.
    xg,wg=np.polynomial.legendre.leggauss(5)
    nodes=(sol.t[1:]+sol.t[:-1])[:,None]/2+np.diff(sol.t)[:,None]*xg/2
    weights=np.diff(sol.t)[:,None]*wg/2
    flow=weights.ravel()@np.array([calc(x,u)[-1] for x,u in zip(nodes.ravel(),sol.sol(nodes.ravel()).T,strict=True)])
    led=abs((sol.y[2,-1]-init[2])-(flow[0]-flow[1]))/init[2]
    assert led<2e-8
    fixed=np.linspace(0.,sol.t[-1],11)
    rows=[]
    for x,u in zip(fixed,sol.sol(fixed).T,strict=True):
        T,Td,actual,peak,rhs,rv,sv,_=calc(x,u)
        rows.append(dict(log_scale_growth=x,T_visible_GeV=T,T_dark_GeV=Td,
                         invariant_light_scalar_yields=u[:2].tolist(),actual_light_scalar_heavy_yields=actual.tolist(),
                         peak_classical_phase_densities=peak.tolist(),net_invariant_rates_over_Hs=rhs[:2].tolist(),
                         dark_to_visible_energy_fraction=u[2]*sv/rv))
    return dict(quadrature=order,stop_max_classical_occupation=closure_limit,include_direct_visible=include_source,
                stopped_at_dilution_threshold=bool(len(sol.t_events[0])),
                energy_input_and_pressure_work_per_visible_entropy=flow.tolist(),
                independent_energy_ledger_relative_defect=led,rows=rows)


def run(eos,seed):
    b=Background(eos)
    row=json.loads(seed.read_text())
    assert abs(row["TF"]-.06)<1e-9 and abs(row["TR"]-.1)<1e-9
    arms=[]
    for use in (False,True):
        controls=[one(b,row,j,.075,use) for j in (0,1)]
        A,B=[m["rows"][-1] for m in controls]
        assert abs(A["T_visible_GeV"]/B["T_visible_GeV"]-1.)<1e-4
        assert np.max(abs(np.array(A["invariant_light_scalar_yields"])/B["invariant_light_scalar_yields"]-1.))<1e-4
        arms.append(dict(direct_visible_source=use,controls=controls))
    return dict(scope="Continuation from an imposed passive seed, not from empty. Selected signed MB dark channels, "
                      "same fast population affinity, visible adiabatic EOS with dark density in H and no feedback "
                      "to the prescribed visible clock. Direct-source case uses FD-only leptonic inflow, no visible "
                      "blocking, dark elastic heating or hadrons. Stops when a phase occupation exceeds .075 as "
                      "a declared limit of the dilute closure, not a physical transition or the final yield.",arms=arms)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("seed",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed),indent=2,allow_nan=False))
