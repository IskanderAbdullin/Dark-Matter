"""Orbit-energy density of states near a finite halo's escape energy.
Static isotropic Vlasov diagnostic only: NOT a Bose collision or evaporation solver.
Requires numpy, scipy and continuation-04/code models.
Run: python continuation-06/code/phase_density_edge.py
"""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent
CODE4=HERE.parents[1] / "continuation-04" / "code"
sys.path.insert(0,str(CODE4))
from critical_lowered import solve as solve_linear
from smooth_lowered_bose import solve2 as solve_quadratic

def phase_density(sol,rt,eta,binding):
    """g(s)=integral_0^r(s) xi^2 sqrt(psi(xi)-s) dxi; no 8*pi² factor."""
    s=float(binding)
    if s>=eta:return 0.
    rmax=rt if s==0 else brentq(lambda r:float(sol.sol(r)[0]-s),1e-5,rt,xtol=3e-13)
    return quad(lambda r:r*r*np.sqrt(max(float(sol.sol(r)[0])-s,0)),
                0,rmax,epsabs=1e-8,epsrel=2e-9,limit=250)[0]

def phase_density_derivative_at_escape(sol,rt):
    # g'(0)=-1/2 integral xi²/sqrt(psi) dxi.
    # xi=rt*(1-u²) regularizes the square-root singularity at the edge.
    def integrand(u):
        if u==0:
            qt=float(sol.sol(rt)[1])
            a=qt/rt**2
            return 2*rt**3/np.sqrt(a*rt)
        r=rt*(1-u*u)
        p=max(float(sol.sol(r)[0]),0)
        return 2*rt*u*r*r/np.sqrt(p)
    return -.5*quad(integrand,0,1,epsabs=1e-8,epsrel=1e-8,limit=200)[0]

def run():
    rows=[]
    for name,solver in (("linear",solve_linear),("quadratic",solve_quadratic)):
        Ic,rt,Q,sol,_=solver(eta=4,nq=128)
        g0=phase_density(sol,rt,4,0)
        gp0=phase_density_derivative_at_escape(sol,rt)
        points=[]
        for s in (0,1e-5,1e-4,.001,.01,.03,.1,.3,.6,1,2,3):
            v=phase_density(sol,rt,4,s)
            points.append(dict(binding=s,g=v,ratio=v/g0))
        for s in (1e-4,.001):
            finite_difference=(phase_density(sol,rt,4,s)-g0)/s
            assert abs(finite_difference/gp0-1)<.01,(name,s)
        assert g0>0 and gp0<0 and points[0]["ratio"]==1
        assert all(points[i]["g"]>points[i+1]["g"] for i in range(len(points)-1))
        rows.append(dict(model=name,eta=4,rt_over_R=rt,Q=Q,Ic=Ic,
                         g0=g0,gp0=gp0,relative_edge_slope=gp0/g0,points=points))
    return dict(models=rows,
      ratio_g0_quadratic_to_linear=rows[1]["g0"]/rows[0]["g0"],
      ratio_mass_quadratic_to_linear=rows[1]["Q"]/rows[0]["Q"],
      checks_passed=8,
      warning="Static state measure only; the microscopic diffusion and true escape flux were not calculated.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
