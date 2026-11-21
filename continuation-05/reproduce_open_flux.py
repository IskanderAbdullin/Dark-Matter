"""Toy one-dimensional open energy-space Bose Fokker-Planck absorbing boundary.
NOT a derived 3D self-gravitating Bose collision operator.
"""
import json
import numpy as np
import math
from scipy.integrate import solve_ivp
from reproduce_cutoff_ladder import occupation

def run(beta=1.,D=1.,flux=.01):
    # s = E_t-E >= 0, J = D*(df/ds-beta*f(1+f)).
    # Fixed J>0 into the escape boundary -> df/ds=beta*f(1+f)+J/D.
    def rhs(s,y):return [beta*y[0]*(1+y[0])+flux/D]
    sol=solve_ivp(rhs,[0,2],[0],rtol=1e-12,atol=1e-14,dense_output=True)
    s=np.array([1e-5,1e-4,1e-3,.01,.1,.5,1,2.])
    f=sol.sol(s)[0]
    j=flux/D; discr=math.sqrt(1-4*j/beta)
    r_plus=(-1+discr)/2; r_minus=(-1-discr)/2
    A0=r_plus/r_minus
    A=A0*np.exp(beta*discr*s)
    f_exact=(r_plus-A*r_minus)/(1-A)
    assert np.max(np.abs(f-f_exact))<3e-11
    pole_s=math.log(1/A0)/(beta*discr)
    assert abs(pole_s-4.679406550517849)<1e-10
    assert abs(f[0]/s[0]/(flux/D)-1)<.0001
    assert abs(f[1]/s[1]/(flux/D)-1)<.0001
    assert np.all(f>0)
    t=np.array([1e-6,1e-5,1e-4])
    edge=[]
    for k in (1,2,3):
        fk=occupation(t,4,k)
        exponent=float(np.polyfit(np.log(t),np.log(fk),1)[0])
        assert abs(exponent-k)<.01
        edge.append(dict(k=k,edge_occupation_power=exponent,
                         inferred_flux_at_edge="finite nonzero" if k==1 else "zero if diffusion D(E_t) finite"))
    return dict(exact_finite_binding_pole_s=pole_s,
      exact_Riccati_solution_max_abs_err=float(np.max(np.abs(f-f_exact))),
      analytical_boundary_slope_fprime_in_s=flux/D,
      rows=[dict(s=float(a),f=float(b)) for a,b in zip(s,f)],
      example_f_at_s_2=float(f[-1]),edge_checks=edge,tests_passed=6,
      important_limitation="Toy 1D Fokker-Planck in energy only. Galactic orbit-averaged coefficients, state densities, Bose kinetic kernels and open-tidal escape not derived.")
if __name__=='__main__':print(json.dumps(run(),indent=2))
