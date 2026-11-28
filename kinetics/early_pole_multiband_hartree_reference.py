"""Multi-band spontaneous/stimulated decay control with shared scalar tadpole."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from early_pole_decay_mode_clock_reference import G,X,KAPPA
from early_pole_shell_hartree_feedback_reference import S,Y,Q


def case(lam,nrad,b=.0,qlow=1e-8):
    ed=np.unique(np.r_[np.geomspace(qlow,30.,nrad+1),1e-5,.1])
    q=np.sqrt(ed[:-1]*ed[1:])
    v=np.diff(ed**3)/(6*np.pi**2)
    vmax=S*Y
    def rhs(t,z,output=False):
        x=X*np.exp(t)
        n=v*np.expm1(np.clip(z,0.,100.))
        nt=n.sum()
        base=1e-9*x
        if lam*nt:
            k=lam/2
            m=brentq(lambda t:t*t-base*base-k*(n@ (1/np.hypot(t,q))),
                     base,base+1.0001*(k*nt)**(1/3)+1e-19,xtol=2e-17,rtol=2e-12)
        else:
            m=base
        if output:
            return dict(x=x,scalar_pole_over_T=m,Y_total=nt/S,
                        p_lt_1e_5_yield=n[ed[1:]<=1e-5].sum()/S,
                        p_lt_point1_yield=n[ed[1:]<=.1].sum()/S,
                        average_p_over_T=float((n@q)/nt) if nt>0. else 0.,
                        peak_occupation=float(np.expm1(np.max(np.clip(z,0.,100.)))),
                        peak_occupation_p_over_T=q[np.argmax(z)],
                        shell_peak_count_p_over_T=q[np.argmax(n/np.diff(np.log(ed)))] if nt>0. else None)
        if m+b>=x:
            return np.zeros_like(q)
        A=(x*x+m*m-b*b)/2
        D=np.sqrt((x*x-(m+b)**2)*(x*x-(m-b)**2))/2
        E=np.hypot(m,q)
        lo=(x*x*E*E+D*D)/(A*E+q*D)
        dx=2*q*D/(m*m)
        hi=lo+dx
        P=D/x
        I=np.log1p(np.exp(-lo)*(-np.expm1(-dx))/(1+np.exp(-hi)))
        ymax=2*P/x*((x+b)**2-m*m)/x**2
        return (KAPPA*x**4*G/(2*E*q*P))*I*ymax
    def ev(t,z):
        return v@np.expm1(np.clip(z,0.,100.))-vmax
    ev.terminal=True;ev.direction=1
    sol=solve_ivp(rhs,(0.,np.log(10.)),np.zeros_like(q),rtol=2e-8,atol=1e-10,
                   method="DOP853",max_step=.035,events=ev,dense_output=True)
    assert sol.success
    checkpoints=[]
    for x in (.0052,.0055,.006,.008,.01,.02,.03,.04,.05):
        t=np.log(x/X)
        if t<=sol.t[-1]+1e-12:
            checkpoints.append(rhs(t,sol.sol(t),True))
    final=rhs(sol.t[-1],sol.y[:,-1],True)
    return dict(lambda_self_quartic=lam,fermion_pole_over_T=b,requested_grid=nrad,momentum_floor=qlow,
                target_hit=bool(len(sol.t_events[0])),points=checkpoints,
                normalized_diagnostics_at_end=final,
                final_n_per_visible_entropy_of_each_bin=list((v*np.expm1(np.clip(sol.y[:,-1],0.,100.)))/S),
                momentum_edges=ed.tolist())


def run():
    result=[case(l,n,b,lo) for l,b in ((0.,0.),(7.65175e-8,0.),
                    (1.08515714e-6,0.),(7.65175e-8,.003))
            for n,lo in ((150,1e-8),(360,1e-9))]
    for i in range(0,len(result),2):
        a,b=result[i]['normalized_diagnostics_at_end'],result[i+1]['normalized_diagnostics_at_end']
        keys=('x','Y_total','p_lt_1e_5_yield','p_lt_point1_yield','average_p_over_T')
        # Avoid claiming convergence in exact arbitrary peak f;
        # in the nonfeedback arm one narrow maximum dominates.
        drift={k:abs(a[k]/b[k]-1.) for k in keys}
        assert max(drift.values())<.045,(i,drift)
        result[i+1]['from_coarser_grid_fractional_drifts']=drift
    return dict(scope="Uniform per-bin occupation and a shared positive scalar tadpole with held FD parent "
                      "and empty fermion daughters across all bins, a prescribed T(a), matched scalar Yukawa "
                      "vacuum width, no complete number/energy closure, collisions or self-consistent "
                      "fermion production. Not a predicted final spectrum or fit.",
                examples=result)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
