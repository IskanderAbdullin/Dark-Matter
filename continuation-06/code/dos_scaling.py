"""Universal eta->0 scaling of the near-escape state density.
Static Bose-Poisson equilibria only; not an evaporation time calculation.
Run: python continuation-06/code/dos_scaling.py
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"continuation-04"/"code"))
from small_eta_universal import density
from critical_lowered import solve
from phase_density_edge import phase_density

def universal():
    eps=1e-5
    def rhs(s,y):return[-y[1]/s**2,s*s*density(y[0])]
    def edge(s,y):return y[0]
    edge.terminal=True;edge.direction=-1
    sol=solve_ivp(rhs,[eps,100],[1-eps**2/6,eps**3/3],
                  events=edge,rtol=1e-10,atol=1e-12,max_step=.03,dense_output=True)
    assert sol.success and len(sol.t_events[0])
    rt=float(sol.t_events[0][0])
    G0=quad(lambda s:s*s*np.sqrt(max(float(sol.sol(s)[0]),0)),
            0,rt,epsrel=1e-10)[0]
    return rt,float(sol.sol(rt)[1]),G0

def run():
    rt,Q,G0=universal()
    rows=[]
    for eta in (.5,.2,.1,.05,.02,.01):
        Ic,r,q,sol,_=solve(eta,nq=168)
        g=phase_density(sol,r,eta,0)
        rows.append(dict(eta=eta,g0=g,g_over_eta2=g/eta**2,
                         relative_error=g/(eta**2*G0)-1))
    assert abs(rows[-1]["relative_error"])<.001
    assert abs(rows[-2]["relative_error"])<.005
    return dict(universal_radius=rt,universal_Q=Q,
       G0=G0,scaling="g(eta*z)~eta**2*G_universal(z)",
       points=rows,checks_passed=2,
       warning="Static state density is not a dynamical escape-rate prediction.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
