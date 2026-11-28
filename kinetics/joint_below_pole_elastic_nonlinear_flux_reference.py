"""Independent nonlinear FD/BE flux checks for the below-pole scalar channel."""
import argparse,json
from pathlib import Path
import numpy as np
from joint_resonant_kinetic_linear_modes_reference import configuration,metric
from joint_below_pole_elastic_linear_modes_reference import setup
from sm_dark_entropy_exchange_reference import StandardModel

def test(r,eos,order,pole):
    T,H,k,N,t,Q,th=configuration(r,eos,6,8.)
    _,L,z=metric(k,th)
    d=z.shape[0]
    eps=1e-5
    p=.14*np.sin(.39*np.arange(d))+.07*np.cos(.31*np.arange(d))
    E,F,S=setup(T,k,t,[x[3] for x in Q],*order,29.,pole,
                perturb=(np.zeros(d),eps*p,-eps*p,p))
    g=(F[1]-F[2])/(2*eps)
    ferr=np.linalg.norm(g-E@p)/np.linalg.norm(E@p)
    err=np.linalg.norm(z.T@F[-1])/np.linalg.norm(F[-1])/np.linalg.norm(z)
    st=abs(S[-1]/(p@F[-1])-1)
    assert np.linalg.norm(F[0])==0 and S[0]==0 and S[-1]>0
    assert ferr<2e-7 and err<1e-12 and st<1e-12,(ferr,err,st)
    return dict(orders=order,pole_cut=pole,total_energy_and_counts_null_error=err,
       independent_small_signal_Jacobian_error=ferr,independent_entropy_projection_error=st,
       nonlinear_small_departure_entropy_per_unit_volume_and_time=S[-1],
       relative_nonlinear_flux_vs_linear=float(np.linalg.norm(F[-1]-E@p)/np.linalg.norm(E@p)))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())["rows"] if r["TR"]==.05 and abs(r["TF"]-.045)<1e-8)
    print(json.dumps(dict(rows=[test(r,eos,o,z) for o,z in (((4,3,2,10),.82),((6,5,3,16),.82),((6,5,3,16),.72))],
       scope="Nonlinear event flux at a manufactured finite-box perturbation; no full physical transport or abundance."),
       indent=2,allow_nan=False))
