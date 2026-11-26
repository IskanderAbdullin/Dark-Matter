"""Particle and heat response of the common scalar-pair collision chart."""

import argparse
import json

import numpy as np

from archived_scalar_pair_reaction_reference import m,mp
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_pair_common_fermion_chart_reference import rate


def extrapolate(output,first):
    large,small=output
    r=(small["finite_difference"]/large["finite_difference"])**2
    J=(np.array(small["particle_heat_response_over_eq_forward_events"])-
       r*np.array(large["particle_heat_response_over_eq_forward_events"]))/(1-r)
    sym=abs(J[0,1]/J[1,0]-1)
    e=np.linalg.eigvalsh((J+J.T)/2)
    assert sym<1e-9 and e[0]>0.
    return dict(richardson_response=J.tolist(),reciprocity_defect=float(sym),
                exact_particle_row_relative_defect=float(np.max(abs(J[0]/np.asarray(first)-1))),
                symmetric_eigenvalues=e.tolist())


def response(T,mu,p,order,steps=(2e-4,6e-5)):
    """Two common-temperature derivatives, in energy units of the anchor T."""
    am=mu/T
    pe=[np.hypot(p,mass) for mass in (mp,m)]
    # x1 is conjugate to particle number and x2 to heat at anchor mu.
    # Log(G/L) = 2*x1 + x2*(E_pair-2*mu)/T for represented thermal events.
    def result(x):
        chem=x[0]-am*x[1]
        bet=np.array([1.-x[1]/2,1.+x[1]/2])/T
        alp=am+np.array([chem,-chem])/2
        th=[np.vstack((r*ee[:-1]-a,r*np.diff(ee))) for r,a,ee in zip(bet,alp,pe,strict=True)]
        S=dict(lab_momentum_edges_GeV=p,entropy_variable_scalar=th[0],entropy_variable_fermion=th[1])
        B,_,F,_=snapshot_functions(S)
        _,_,gross,closure,_=rate(p,F,B,order)
        diff=gross[0]-gross[1]
        y=np.array([2*diff[0],diff[1]/T-2*am*diff[0]])
        return y,gross,closure

    equilibrium,w,cons=result(np.zeros(2))
    reference=w[1,0]
    tiny=np.max(np.abs(equilibrium))/reference
    assert tiny<1e-12
    predicted_first=np.array([4*reference,2*(w[1,1]/T-2*am*reference)])/reference
    output=[]
    for h in steps:
        col=[]
        for dx in np.eye(2)*h:
            yp,_,_=result(dx)
            ym,_,_=result(-dx)
            col.append((yp-ym)/(2*h*reference))
        J=np.array(col).T
        err=np.max(abs(J[0]-predicted_first)/predicted_first)
        sym=abs(J[0,1]/J[1,0]-1)
        eig=np.linalg.eigvalsh((J+J.T)/2)
        assert err<2e-7 and sym<2e-7 and min(eig)>0.
        output.append(dict(finite_difference=h,particle_heat_response_over_eq_forward_events=J.tolist(),
                           theoretical_first_row_error=float(err),reciprocity_defect=float(sym),
                           eigenvalues_of_symmetric_part=eig.tolist()))
    return dict(T_GeV=T,mu_GeV=mu,mesh_cells=len(p)-1,quadrature=order,
                event_check=cons,initial_particle_heat_residual_per_forward_event=float(tiny),
                first_row_from_analytic_event_moments=predicted_first.tolist(),
                equilibrium_forward_event_number_energy=w[1].tolist(),controls=output,
                extrapolated=extrapolate(output,predicted_first))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fine",action="store_true")
    a=p.parse_args()
    order=(4,10,10,4) if a.fine else (3,7,6,2)
    grid=np.r_[0.,np.geomspace(.00004,.0015,6),np.geomspace(.003,1.07,11)]
    print(json.dumps(response(.017,-.030,grid,order),indent=2,allow_nan=False))
