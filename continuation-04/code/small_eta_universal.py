"""eta->0 universal shallow-cutoff exactly critical lowered Bose equilibrium.
Requires numpy scipy and critical_lowered.py in same directory.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from critical_lowered import solve
nodes,weights=leggauss(128)
t=(nodes+1)/2;w=weights/2
def density(U):
    U=float(np.clip(U,0,1))
    return 0 if U==0 else 1.5*U**2.5*np.dot(w,t*t*(1-t*t)/(1-U*(1-t*t)))
def run():
    eps=1e-5
    def ode(s,a):return[-a[1]/s**2,s*s*density(a[0])]
    def edge(s,a):return a[0]
    edge.terminal=True;edge.direction=-1
    sol=solve_ivp(ode,[eps,100],[1-eps**2/6,eps**3/3],
        events=edge,rtol=1e-10,atol=1e-12,max_step=.03,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1
    st=float(sol.t_events[0][0]);Qt=float(sol.sol(st)[1])
    rows=[]
    for eta in (.5,.2,.1,.05,.02,.01):
        Ic,rt,Q,_,_=solve(eta,nq=168)
        rows.append(dict(eta=eta,rt_over_sqrt_eta=rt/np.sqrt(eta),
            Qt_over_eta_3_over_2=Q/eta**1.5,
            slope_times_sqrt_eta=np.pi*np.sqrt(eta)/(2*np.sqrt(6)*Ic)))
    slope=3*np.pi/(4*np.sqrt(6))
    assert abs(rows[-1]['rt_over_sqrt_eta']/st-1)<.001
    assert abs(rows[-1]['Qt_over_eta_3_over_2']/Qt-1)<.001
    assert abs(rows[-1]['slope_times_sqrt_eta']/slope-1)<.001
    return dict(universal_rt=st,universal_Qt=Qt,slope_limit=slope,
                finite_eta=rows,checks_passed=3,
                caveat="Static mathematical limit, no formation or stability theorem.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
