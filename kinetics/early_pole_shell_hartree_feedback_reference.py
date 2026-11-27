"""Self-induced scalar tadpole and a restricted early shell occupation clock."""
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import roots_legendre
from early_pole_decay_mode_clock_reference import kernel, X, Q, KAPPA

S=2*np.pi**2/45*106.75
Y=.005056924750678494
ZGOAL=S*Y*.92


def setup(qs=(.5,1.5),N=10):
    x,w=roots_legendre(N)
    lo,hi=np.array(qs)*Q
    p=(lo+hi)/2+(hi-lo)/2*x
    ww=(hi-lo)/2*w*p*p/(2*np.pi**2)
    vol=(hi**3-lo**3)/(6*np.pi**2)
    assert abs(ww.sum()/vol-1)<1e-14
    return p,ww/vol,vol


def pole(n,lam,base,p,w):
    if lam==0 or n==0:
        return base
    pref=lam*n/2
    top=base+1.001*pref**(1/3)+1e-20
    z=brentq(lambda m:(m*m-base*base)-pref*(w@ (1/np.hypot(p,m))),
             base,top,xtol=2e-17,rtol=2e-12)
    return z


def case(lam,b,mu,qbounds=(.5,1.5),rad=10):
    p,w,phase=setup(qbounds,rad)
    ugoal=np.log1p(ZGOAL/phase)
    def coeff(x,u):
        # Stop before evaluating any trial that represents an unphysical
        # arbitrarily large occupation in the artificial infinite parent reservoir.
        n=phase*np.expm1(np.clip(u,0.,ugoal+4.))
        m=pole(n,lam,1e-9*x,p,w)
        if x<=m+b:
            return 0.,0.,m,n
        k=[kernel(x,m,b,float(q),mu) for q in p]
        aa=w@np.array([q['seed_over_width']*q['vacuum_matched_yukawa_width_ratio'] for q in k])
        bb=w@np.array([q['stimulus_over_width']*q['vacuum_matched_yukawa_width_ratio'] for q in k])
        return KAPPA*x*x*aa,KAPPA*x*x*bb,m,n
    def f(t,u):
        a,b0,_,_=coeff(X*np.exp(t),u[0])
        return [a*np.exp(-max(u[0],0.))+b0*(-np.expm1(-max(u[0],0.)))]
    def reached(t,u):
        return u[0]-ugoal
    reached.terminal=True;reached.direction=1
    sol=solve_ivp(f,(0.,np.log(10.)),[0.],rtol=5e-9,atol=1e-11,
                  max_step=.018,dense_output=True,events=reached)
    assert sol.success
    rows=[]
    for x in (.005,.0055,.006,.008,.01,.02,.05):
        if np.log(x/X)>sol.t[-1]+1e-10:
            continue
        a,b0,m,n=coeff(x,sol.sol(np.log(x/X))[0])
        rows.append(dict(parent_over_T=x,scalar_over_T=m,shell_n_over_T3=n,
                         shell_yield=n/S,seed_over_H=a,stimulus_over_H=b0))
    a,b0,m,n=coeff(X*np.exp(sol.t[-1]),sol.y[0,-1])
    return dict(lambda_self_quartic=lam,fermion_pole_over_T=b,held_fermion_mu_over_T=mu,
                shell_edges_relative_to_pivot=qbounds,completed_to_x=float(X*np.exp(sol.t[-1])),
                reached_nominal_cold_count=bool(len(sol.t_events[0])),
                final_scalar_pole_over_T=m,final_shell_yield=n/S,
                final_count_over_goal=n/ZGOAL,final_seed_over_H=a,points=rows)


def run():
    low=7.65175e-8
    high=1.08515714e-6
    # Equilibrium check for the conventional vacuum-subtracted, one-real-scalar
    # massless tadpole. This is NOT what is assigned to the initially empty shell.
    eq=quad(lambda q:q*np.exp(-q)/(-np.expm1(-q)),0.,100.)[0]/(2*np.pi**2)
    assert abs(eq*12-1)<1e-12
    p,w,_=setup()
    checks=[]
    for lam in (low,high):
        for n in (1e-13,1e-8,1e-4,ZGOAL):
            m=pole(n,lam,1e-9*X,p,w)
            defect=abs((m*m-(1e-9*X)**2)/(lam*n/2*(w@(1/np.hypot(p,m))))-1)
            assert defect<7.e-8
            checks.append(defect)
    rows=[case(l,b,mu) for l,b,mu in ((0.,0.,None),(low,0.,None),(high,0.,None),
            (low,.003,None),(high,.003,None),(low,0.,-3e-5))]
    k=case(low,0.,None,rad=18)
    assert abs(k['final_shell_yield']/rows[1]['final_shell_yield']-1)<.02
    alt=case(low,0.,None,qbounds=(.8,1.2))
    return dict(scope="One top-hat comoving scalar shell, imposed radiation expansion T proportional to 1/a "
                      "and fixed parent FD pole; no energy closure, no evolving fermion daughter, no other "
                      "scalar momenta or spectral-width effects. Only illustrative Hartree tadpole, "
                      "not a real early-Universe self-energy or bound on the full production network.",
                shell_nominal_target_Y=.92*Y,independent_thermal_tadpole_integral=eq,
                max_gap_equation_residual=max(checks),examples=rows,
                independent_refined_quadrature_case=k,alternative_narrow_shell_case=alt)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
