"""Conditional nonlinear e/mu/tau vacuum MB portal in two thermal baths."""
import argparse
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from mb_number_heat_state_reference import state
from sm_dark_entropy_exchange_reference import StandardModel, MPL
from timelike_chemical_heat_eigenmodes_reference import rates, trial

T0=.6882487043

class Portal:
    def __init__(self, nodes, moments):
        beta=1/nodes[::-1]
        nu,q,q2=moments[::-1].T
        self.lo,self.hi=beta[[0,-1]]
        self.f=CubicHermiteSpline(beta,np.log(nu),-q/nu,extrapolate=False)
        # Curvature is linear on each interval, so checking the actual
        # two one-sided endpoint values suffices for convexity everywhere.
        a,b=self.f.c[:2]
        self.min_curv=float(min(np.min(2*b),np.min(6*a*np.diff(beta)+2*b)))
        assert self.min_curv>0 and np.all(q2/nu-(q/nu)**2>0)
        self.mean=self.f.derivative()

    def flow(self,Ts,Td,eta):
        s,d=1/Ts,1/Td
        assert self.lo<=min(s,d) and max(s,d)<=self.hi
        ls,ld=self.f(s),self.f(d)+2*eta
        es,ed=-self.mean(s),-self.mean(d)
        # Subtract positive event fluxes without losing the linear response
        # at equal temperatures or multiplying an overflowing fugacity.
        if ld>=ls:
            ref=np.exp(ld)
            r=-ref*np.expm1(ls-ld)
            q=ed*r+np.exp(ls)*(ed-es)
        else:
            ref=np.exp(ls)
            r=ref*np.expm1(ld-ls)
            q=es*r+np.exp(ld)*(ed-es)
        return r,q,2*eta*r+(s-d)*q

def equilibrium_n(T):
    return np.exp(state(T,T**3)["log_eq_density"])

def run_case(eos,portal,dark_ratio,number_ratio,*,scale=1.,rtol=5e-11,form="temperature"):
    Td0=T0*dark_ratio
    n0=equilibrium_n(Td0)*number_ratio
    sm0=eos.state(T0)
    dark0=state(Td0,n0)
    init=np.r_[[0.,0.,0.],np.zeros(6)]
    def physical(x,z):
        n=n0*np.exp(z[2]-3*x)
        if form=="temperature":
            return T0*np.exp(z[0]),Td0*np.exp(z[1]),n
        s=sm0[2]*np.exp(z[0]-3*x)
        sd=dark0["entropy"]*np.exp(z[1]-3*x)
        Ts=brentq(lambda t:eos.state(t)[2]/s-1,np.exp(eos.lo)*1.000000000001,
                  np.exp(eos.hi)*.999999999999,xtol=1e-14)
        Td=brentq(lambda t:state(t,n)["entropy"]/sd-1,.002,1.4,xtol=1e-14)
        return Ts,Td,n
    def end(x,z):
        return physical(x,z)[0]-.01
    end.terminal=True
    end.direction=-1

    def rhs(x,z):
        Ts,Td,n=physical(x,z)
        rs,ps,ss,cs=eos.state(Ts)
        dark=state(Td,n)
        rd,pd,sd,cv,e,eta=(dark[k] for k in ("rho","pressure","entropy",
                           "cv_fixed_number","chemical_energy","eta"))
        H=np.sqrt((rs+rd)/3)/MPL
        R,Q,P=scale/H*np.array(portal.flow(Ts,Td,eta))
        assert P>=-3e-22
        a3,a4=np.exp(3*x),np.exp(4*x)
        if form=="temperature":
            ds=(-3*(rs+ps)+Q)/(cs*Ts)
            dd=(-3*(rd+pd-n*e)-Q+2*e*R)/(cv*Td)
        else:
            ds=Q/(Ts*ss)
            dd=(-Q/Td+2*eta*R)/sd
        return [ds,dd,-2*R/n,a4*(rs-3*ps),a4*(rd-3*pd),a4*Q,a3*P,
                -2*a3*R,a3*Q/Ts]

    sol=solve_ivp(rhs,(0.,8.),init,method="DOP853",rtol=rtol,atol=rtol*.003,
                  max_step=.009,events=end,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1
    x=sol.t[-1]
    z=sol.y[:,-1]
    Ts,Td,n=physical(x,z)
    rs,ps,ss,_=eos.state(Ts)
    dark=state(Td,n)
    a3,a4=np.exp(3*x),np.exp(4*x)
    checks=dict(energy=max(abs(a4*rs-sm0[0]-z[3]-z[5]),
                           abs(a4*dark["rho"]-dark0["rho"]-z[4]+z[5]))/(sm0[0]+dark0["rho"]),
                entropy=abs(a3*(ss+dark["entropy"])-sm0[2]-dark0["entropy"]-z[6])/(sm0[2]+dark0["entropy"]),
                number=abs(a3*n-n0-z[7])/n0,
                sm_entropy=abs(a3*ss-sm0[2]-z[8])/sm0[2])
    assert max(checks.values())<1e-6, (form,dark_ratio,number_ratio,scale,checks)
    assert z[6]>=0
    xx=np.linspace(0,x,101)
    samples=[]
    flow=[]
    for u,y in zip(xx,sol.sol(xx).T):
        t,v,w=physical(u,y)
        eta=state(v,w)["eta"]
        samples.append((t,v/t,w*np.exp(3*u)/n0,eta))
        flow.append(portal.flow(t,v,eta)[:2])
    samples=np.asarray(samples)
    flow=np.asarray(flow)*scale
    if scale==0.:
        assert abs(a3*n/n0-1)<2e-12 and abs(a3*ss/sm0[2]-1)<5e-7
        assert abs(a3*dark["entropy"]/dark0["entropy"]-1)<5e-7
    return dict(form=form,scale=scale,initial_dark_to_SM_ratio=dark_ratio,
                initial_over_equilibrium_number=number_ratio,final_T_ratio=Td/Ts,
                final_comoving_number_fraction=a3*n/n0,final_eta=dark["eta"],
                final_log_a=x,relative_ledger_errors=checks,entropy_produced=z[6],
                number_flux_positive_samples=int(np.sum(flow[:,0]>0)),
                heat_flux_positive_samples=int(np.sum(flow[:,1]>0)),nfev=sol.nfev,
                min_max_sampled_T_ratio=[float(np.min(samples[:,1])),float(np.max(samples[:,1]))])

def test_flows(p, eos):
    pairs=[(.05,.05),(.3,.299999),(.29,.31),(.688,.69),(.6,.8),(.8,.6)]
    disorder=[]
    for Ts,Td in pairs:
        for eta in (-.4,-.001,0.,.001,.4):
            r,q,prod=p.flow(Ts,Td,eta)
            r2,q2,prod2=p.flow(Td,Ts,-eta)
            assert prod>=-1e-35
            # Reverse swaps the baths but the rate normalization changes
            # because only the first bath in flow is fugacity weighted.
            ref=np.exp(2*eta)
            for a,b in ((r,-ref*r2),(q,-ref*q2),(prod,ref*prod2)):
                err=abs(a-b)/max(abs(a),abs(b),1e-200)
                assert err<5e-8
                disorder.append(err)
    from timelike_chemical_heat_eigenmodes_reference import thermo
    T=.6882487043
    sm=eos.state(T)
    chi,rd=thermo(T)
    H=np.sqrt((sm[0]+rd)/3)/MPL
    sd=np.linalg.inv(chi)
    def f(y):
        eta,dt=sd@y
        Ts=T*(1-y[1]/sm[3])
        Td=T*(1+dt)
        r,q,_=p.flow(Ts,Td,eta)
        return np.array([2*r,q/T])
    base=chi[0,0]*2e-6
    J=np.stack([(f(v)-f(-v))/(2*base*H) for v in base*np.eye(2)],axis=1)
    ev=np.sort(np.linalg.eigvals(J).real)
    target=np.array(trial(T,eos)["eigenmodes_per_H"])
    err=np.max(abs(ev/target-1))
    assert err<.003, (ev,target,err)
    return dict(max_bath_reversal_relative_error=max(disorder),
                linear_modes_per_H=ev.tolist(),linear_modes_reference_per_H=target.tolist(),
                linear_modes_max_relative_error=err)


def main(table, grid):
    eos=StandardModel(table)
    data=np.load(grid)
    t,m=data[:,0],data[:,1:]
    p=Portal(t,m)
    coarse=Portal(t[::2],m[::2])
    errs=[]
    for i in range(3,len(t)-3,8):
        b=1/t[i]
        y=m[i]
        errs.append(np.max(abs(np.exp(coarse.f(b))*np.r_[1.,-coarse.mean(b)]/y[:2]-1)))
    assert p.flow(T0,T0,0.)==(0.,0.,0.)
    cases=[(1.,1.),(.75,.6),(1.25,1.6),(.85,1.25),(1.15,.7)]
    rows=[run_case(eos,p,*pair) for pair in cases]
    zero=[run_case(eos,p,*pair,scale=0.) for pair in cases[:3]]
    comparison=run_case(eos,p,*cases[0],form="entropy",rtol=1e-11)
    reduced=run_case(eos,coarse,*cases[0],rtol=1e-11)
    effects={}
    for k in ("final_T_ratio","final_comoving_number_fraction"):
        effects[k]=dict(temperature_vs_entropy=abs(rows[0][k]-comparison[k]),
                        refined_vs_coarse_rate=abs(rows[0][k]-reduced[k]))
        assert max(effects[k].values())<3e-4
    return dict(rows=rows,zero=zero,entropy_control=comparison,coarse_control=reduced,
                exact_symmetry_and_linear_validation=test_flows(p,eos),
                max_midpoint_rate_moment_relative_error=max(errs),grid_points=len(t),
                minimum_log_rate_convexity=p.min_curv,convergence=effects,
                scope="Conditional symmetric MB dark partner plus bath scalar, SM EOS, direct e/mu/tau timelike vacuum only; not a physical relic calculation.")

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("eos")
    a.add_argument("grid")
    a.add_argument("--build-grid",action="store_true")
    a.add_argument("--nodes",type=int,default=161)
    args=a.parse_args()
    if args.build_grid:
        ts=np.geomspace(.0025,1.2,args.nodes)
        np.save(args.grid,np.c_[ts,[rates(t) for t in ts]])
    else:
        print(json.dumps(main(args.eos,args.grid),indent=2,allow_nan=False))
