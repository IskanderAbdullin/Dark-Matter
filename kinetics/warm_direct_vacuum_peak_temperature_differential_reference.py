"""Differentiate yields with respect to the upper endpoint, not a reheating dynamics fit."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from kloe_pion_band_portal_reference import table
from sm_dark_entropy_exchange_reference import MPL
from warm_late_cold_distance_source_reference import q_limit
from warm_late_cold_visible_exact_recoil_fold_reference import local
from warm_late_cold_visible_kloe_recoil_reference import rates
from warm_light_occupation_stage_audit_reference import read


def run(eos,kloe,clock):
    bg=Background(eos)
    data=table(kloe)
    C=read(clock)
    x,w=roots_legendre(8)
    out=[]
    T0=.05
    for i,D in enumerate((None,.3,1.)):
        y=sum(C["cases"][i]["controls"][-1]["standard_FD_lepton_and_BE_measured_pion_yield"])
        def source(t):
            z=np.log(t)
            a=bg.scale(t)
            q=.1*bg.T0 if D is None else q_limit(bg,a,D)[0]
            pc=q/a
            h=t*t*np.sqrt(np.pi**2/90*np.exp(bg.gr(z)))/MPL
            s=2*np.pi**2/45*t**3*np.exp(bg.gs(z))
            d=(3+bg.gs.derivative()(z))/3
            return (local(t,pc,12)[0]+rates(t,pc,*data,12)[1])*d/(s*h)
        exact=source(T0)/y
        controls=[]
        for spread in (2e-3,1e-3,5e-4):
            lo=np.log(T0)-spread;hi=np.log(T0)+spread
            direct=spread*(w@np.asarray([source(t) for t in np.exp((hi+lo)/2+spread*x)]))/(2*spread)
            controls.append(dict(log_T_half_step=spread,dlnY_dln_peak=float(direct/y)))
        assert abs(controls[-1]["dlnY_dln_peak"]/exact-1)<.001
        diffs=[]
        for dT in (-.002,-.001,-.0005,.0005,.001,.002):
            lo,hi=np.log([T0,T0+dT])
            du=(hi-lo)/2
            dY=du*(w@np.array([source(t) for t in np.exp((hi+lo)/2+du*x)]))
            diffs.append(dict(delta_peak_MeV=dT*1000,yield_ratio_to_peak50=float(1+dY/y),
                              first_order_ratio_at_peak=float(1+exact*dT/T0)))
        out.append(dict(selection=C["cases"][i]["selection"],initial_yield=y,
                        analytic_log_endpoint_derivative=float(exact),centered_temperature_checks=controls,
                        finite_peak_shifts=diffs))
    return dict(scope="Different upper temperature cut in the same prescribed adiabatic radiation source, "
                      "zero other changes; not a response to an actual reheating perturbation or solved alternate expansion.",
                cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","clock"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.clock),indent=2,allow_nan=False))
