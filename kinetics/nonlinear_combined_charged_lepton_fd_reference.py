"""Combined conservative FD benchmark with separate spacelike and timelike kernels."""
import argparse
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator,PchipInterpolator
from boost_axis_elastic_heat_reference import integrate as boost_integrate, leptons
from nonlinear_timelike_fd_trajectory_reference import run
from positive_fd_portal_quadrature_reference import PositivePortal
from sm_dark_entropy_exchange_reference import StandardModel

TS=np.array([.0094,.015,.027,.05,.085,.13,.2,.29,.42,.59,.77,1.])
RATIOS=np.array([.37,.49,.62,.77,.9,1.,1.12,1.25,1.45])
ETAS=np.array([-1.7,-.8,-.35,0.,.35,.75,1.3,2.4,4.1,6.2])

def coef(T,r,eta,order=14,cut=28):
    q=sum(boost_integrate(T,T*r,eta,m,n=order,cut=cut,linear=r==1.) for m in leptons)
    return q/T**9 if r==1. else q/(1-r)/T**9

def build(file,refine=False):
    T=np.sort(np.r_[TS,.105,.168,.24,.35,.5,.67]) if refine else TS
    ratios=np.sort(np.r_[RATIOS,.55,.705,.835,.95,1.06,1.19]) if refine else RATIOS
    etas=np.sort(np.r_[ETAS,-.575,-.175,.175,.55,1.025,1.85]) if refine else ETAS
    values=np.array([[[coef(t,r,z,12 if not refine else 14) for z in etas]
                      for r in ratios] for t in T])
    assert values.min()>0 and np.all(np.isfinite(values))
    np.savez(file,T=T,r=ratios,z=etas,v=values)

class Combined:
    def __init__(self,file,linear_T=False):
        d=np.load(file)
        self.T,self.r,self.eta=np.log(d["T"]),d["r"],d["z"]
        v=np.log(d["v"])
        self.ip=PchipInterpolator(self.T,v,axis=0,extrapolate=False)
        self.linear_T=linear_T
        self.alt=RegularGridInterpolator((self.T,self.r,self.eta),v)
        self.p=PositivePortal()
    def coefficient(self,T,r,z):
        assert self.T[0]<=np.log(T)<=self.T[-1] and self.r[0]<=r<=self.r[-1] and self.eta[0]<=z<=self.eta[-1],(T,r,z)
        if self.linear_T:
            return np.exp(self.alt([np.log(T),r,z])[0])
        i=min(self.r.size-2,np.searchsorted(self.r,r,side="right")-1)
        j=min(self.eta.size-2,np.searchsorted(self.eta,z,side="right")-1)
        u=(r-self.r[i])/(self.r[i+1]-self.r[i])
        v=(z-self.eta[j])/(self.eta[j+1]-self.eta[j])
        return np.exp(np.array([1-u,u])@self.ip(np.log(T))[i:i+2,j:j+2]@np.array([1-v,v]))
    def flow(self,Ts,Td,eta):
        z=self.p.flow(Ts,Td,eta)
        r=Td/Ts
        # Run expects Q positive from dark to SM. Both operators
        # use the same sign convention in this sum.
        q=(r-1)*Ts**9*self.coefficient(Ts,r,eta)
        return z+np.array([0.,q,q*(Td-Ts)/(Td*Ts)])

def main(eosfile,grid,quick=False):
    p=Combined(grid)
    eos=StandardModel(eosfile)
    errors=[]
    alternate=Combined(grid,linear_T=True)
    old_errors=[]
    for t,r,z in ((.013,.52,2.4),(.05,.56,.8),(.1,.65,.2),
                  (.22,.7,.45),(.36,.88,.2),(.52,.95,.03),
                  (.64,.985,.025),(.72,1.15,.3)):
        c=coef(t,r,z,22,34)
        b=p.coefficient(t,r,z)
        err=abs(c/b-1)
        assert err<.11,(t,r,z,err)
        errors.append(err)
        old_errors.append(abs(c/alternate.coefficient(t,r,z)-1))
    local=run(eos,p)
    comp=run(eos,PositivePortal())
    ctrl=None if quick else run(eos,alternate,rtol=7e-9)
    shifts={k:local[k]-comp[k] for k in
            ("T_ratio_final","comoving_N_final_over_initial","final_n_over_sm_entropy")}
    for t in (.05,.3,.7):
        assert np.all(p.flow(t,t,0.)==0)
    d=None if ctrl is None else {k:abs(local[k]-ctrl[k]) for k in shifts}
    if d is not None: assert max(d.values())<.15*max(abs(k) for k in shifts.values()),d
    return dict(combined=local,timelike_only=comp,elastic_increment=shifts,
               linear_temperature_control=ctrl,interpolation_path_differences=d,
               maximum_linear_T_interp_fraction=max(old_errors),
               maximum_elastic_interp_fraction=max(errors),
               scope="Conditional charged-lepton-only FD effective vacuum portal; no medium corrections, real vector, continuous hadrons, or physical final relic.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("grid")
    p.add_argument("--build",action="store_true")
    p.add_argument("--refine",action="store_true")
    p.add_argument("--quick",action="store_true")
    a=p.parse_args()
    if a.build: build(a.grid,a.refine)
    else: print(json.dumps(main(a.eos,a.grid,a.quick),indent=2,allow_nan=False))
