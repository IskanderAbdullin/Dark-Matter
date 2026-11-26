"""Two kinetic temperatures and evolving symmetric FD pair density.

Restricted benchmark: finite-width vacuum e/mu/tau portal; no spacelike
heat, real vector particles, thermal propagators or full hadronic input.
"""
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from fd_number_heat_state_reference import state,evaluate
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from nonlinear_timelike_mb_trajectory_reference import (
    Portal as MBPortal, run_case as mb_case, equilibrium_n, T0
)
from positive_fd_portal_quadrature_reference import PositivePortal

def run(eos,portal,ratio=1.,eta=0.,*,form="temperature",rtol=1e-8,scale=1.):
    Td0=ratio*T0
    init=evaluate(Td0,eta)
    SM=eos.state(T0)
    n0=init["n"]
    rnorm=SM[0]+init["rho"]
    snorm=SM[2]+init["entropy"]
    def physical(x,y):
        n=n0*np.exp(y[2]-3*x)
        if form=="temperature":
            return T0*np.exp(y[0]),Td0*np.exp(y[1]),n
        s=SM[2]*np.exp(y[0]-3*x)
        d=init["entropy"]*np.exp(y[1]-3*x)
        Ts=brentq(lambda t:eos.state(t)[2]/s-1,np.exp(eos.lo)*1.0000000001,
                   np.exp(eos.hi)*.9999999999,xtol=1e-14)
        Td=brentq(lambda t:state(t,n)["entropy"]/d-1,.002,1.35,xtol=1e-14)
        return Ts,Td,n
    def stop(x,y):
        return physical(x,y)[0]-.01
    stop.terminal=True
    stop.direction=-1

    def rhs(x,y):
        Ts,Td,n=physical(x,y)
        rs,ps,ss,cs=eos.state(Ts)
        d=state(Td,n)
        rd,pd,sd,cv,en,eta_d=(d[k] for k in ("rho","pressure","entropy",
                                 "cv_fixed_number","chemical_energy","eta"))
        H=np.sqrt((rs+rd)/3)/MPL
        R,Q,P=scale/H*np.array(portal.flow(Ts,Td,eta_d)) if scale else (0.,0.,0.)
        assert P>=-1e-200
        if form=="temperature":
            ssx=(-3*(rs+ps)+Q)/(cs*Ts)
            sdx=(-3*(rd+pd-n*en)-Q+2*en*R)/(cv*Td)
        else:
            ssx=Q/(Ts*ss)
            sdx=(2*eta_d*R-Q/Td)/sd
        a3,a4=np.exp(3*x),np.exp(4*x)
        return [ssx,sdx,-2*R/n,a4*(rs-3*ps)/rnorm,a4*(rd-3*pd)/rnorm,
                a4*Q/rnorm,a3*P/snorm,-2*a3*R/n0,a3*Q/(Ts*snorm)]
    runs=solve_ivp(rhs,(0.,7.),np.zeros(9),method="Radau",rtol=rtol,atol=rtol*.02,
                   max_step=.045,events=stop)
    assert runs.success and len(runs.t_events[0])==1,runs.message
    x=runs.t[-1]
    y=runs.y[:,-1]
    Ts,Td,n=physical(x,y)
    rs,ps,ss,_=eos.state(Ts)
    d=state(Td,n)
    a3,a4=np.exp(3*x),np.exp(4*x)
    checks=dict(SM_energy=abs((a4*rs-SM[0])/rnorm-y[3]-y[5]),
                dark_energy=abs((a4*d["rho"]-init["rho"])/rnorm-y[4]+y[5]),
                total_entropy=abs(a3*(ss+d["entropy"])/snorm-1-y[6]),
                SM_entropy=abs((a3*ss-SM[2])/snorm-y[8]),
                number=abs(a3*n/n0-1-y[7]))
    assert max(checks.values())<2e-6 and y[6]>=0.,checks
    if scale==0.:
        assert abs(a3*n/n0-1)<1e-13 and abs(a3*d["entropy"]/init["entropy"]-1)<2e-6
    return dict(T_initial=T0,T_ratio_initial=ratio,eta_initial=eta,n_initial=n0,
            T_ratio_final=Td/Ts,T_SM_final=Ts,eta_final=d["eta"],comoving_N_final_over_initial=a3*n/n0,
            final_n_over_sm_entropy=n/ss,entropy_production=y[6]*snorm,log_a_final=x,
            ledger_residuals=checks,coordinate_form=form,nfev=runs.nfev,scale=scale)

def main(eosfile,mbgrid,quick=False):
    eos=StandardModel(eosfile)
    p=PositivePortal()
    runs=[]
    for ratio,eta in ((1.,0.),) if quick else ((1.,0.),(.8,-.4),(1.2,.35)):
        runs.append(run(eos,p,ratio,eta))
    base=runs[0]
    no=run(eos,p,scale=0.)
    out=dict(runs=runs,zero_coupling=no,
        scope="Nonlinear conditional FD e/mu/tau vacuum-only trajectory; no full relic or complete hadronic input.")
    data=np.load(mbgrid)
    pmb=MBPortal(data[:,0],data[:,1:])
    # Match initial temperature and actual total particle density, not
    # separately the different equilibrium densities of FD and MB.
    z=mb_case(eos,pmb,1.,base["n_initial"]/equilibrium_n(T0))
    out["MB_control_at_same_initial_number"]=z
    if not quick:
        en=run(eos,p,form="entropy",rtol=3e-9)
        fine=run(eos,PositivePortal((12,34,22)),rtol=3e-9)
        diff={k:[abs(base[k]-v[k]) for v in (en,fine)] for k in
                 ("T_ratio_final","comoving_N_final_over_initial","final_n_over_sm_entropy")}
        assert max(max(x) for x in diff.values())<1e-4,diff
        out.update(entropy_crosscheck=en,fine_quadrature=fine,differences=diff)
    return out

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("eos")
    a.add_argument("mbgrid")
    a.add_argument("--quick",action="store_true")
    args=a.parse_args()
    print(json.dumps(main(args.eos,args.mbgrid,args.quick),indent=2,allow_nan=False))
