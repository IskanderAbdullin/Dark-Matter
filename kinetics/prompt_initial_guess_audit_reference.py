"""Check that a different nonlinear guess leaves the implicit target unchanged."""
import argparse,json
from pathlib import Path
import numpy as np
from joint_resonant_centered_stiff_reference import Centered,advance
from joint_resonant_kinetic_linear_modes_reference import configuration,metric
from sm_dark_entropy_exchange_reference import StandardModel


def test(eos,r):
    T,H,k,N,t,Q,th=configuration(r,eos,8,6.5)
    x=np.concatenate([a.T.ravel() for a in th])
    x+=.03*np.cos(np.arange(len(x))*.44)
    op=Centered(k,H)
    U0,_s,_c=op.moment(x)
    h=.01
    d,a=advance(op,x,h,T,vacuum_scaled=True)
    logs=[]
    q=None
    diagnostic={}
    for qh in (h/64,h/16,h/4,h):
        q,b=advance(op,x,qh,T,vacuum_scaled=True,initial_theta=q,
                    observer=lambda i,r,j:logs.append((i,r,j)),centered_diagnostic=diagnostic)
    Ua=op.moment(d)[0]
    Ub=op.moment(q)[0]
    err=float(np.max(abs(Ua-Ub)/np.repeat(U0[::2],2)))
    z0,rat,_=op.reference(x,T)
    f=op.centered(diagnostic["centered"],diagnostic["reference"],diagnostic["parent_ratio"])[0]
    residual=float(np.max(abs(Ub-U0-h*f)/np.repeat(np.maximum(U0[::2],Ub[::2]),2)))
    f_round=op.centered(q-z0,z0,rat)[0]
    rounded=float(np.max(abs(Ub-U0-h*f_round)/np.repeat(np.maximum(U0[::2],Ub[::2]),2)))
    _,_,inv=metric(k,th)
    invariance=float(np.max(abs((inv.T@Ub)/(inv.T@U0)-1.)))
    assert err<2e-10 and residual<8.1e-12 and invariance<1e-10 and len(logs)>=4,(err,residual,invariance,len(logs),b["residual"])
    return dict(ordinary_final_residual=a["residual"],continued_same_target_residual=residual,
                ordinary_vs_continued_moment_difference=err,invariant_defect=invariance,
                observed_solver_iteration_count=len(logs),rounded_theta_only_residual=rounded)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args()
    r=next(x for x in json.loads(Path(a.spectra).read_text())["rows"] if x["TR"]==.05 and abs(x["TF"]-.045)<1e-8)
    print(json.dumps(test(StandardModel(a.eos),r),indent=2,allow_nan=False))
