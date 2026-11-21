#!/usr/bin/env python3
"""Independent equilibrium Bose-halo Poisson solver and sanity tests.

This does NOT refit rotation-curve observations or assess cosmological viability.
Source equation: Abdullin & Popov, JCAP11(2021)055.
Requires Python 3.10+, numpy, scipy, mpmath; matplotlib for graphs.
"""
from pathlib import Path
import json
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.constants import G, hbar, c, e as elementary_charge, parsec

mp.mp.dps = 24
NU_CRIT = float(mp.zeta(1.5))

def central_fugacity(nu):
    assert 0 < nu <= NU_CRIT
    if abs(nu - NU_CRIT) < 1.e-12:
        return 1.0
    return brentq(lambda z: float(mp.polylog(1.5,z))-nu,
                  min(1.e-15,nu/10), 1-1.e-15, xtol=1.e-15)

def equilibrium(nu, alpha=0.0, gamma=1.0, xi_max=12.0):
    """x_n=Li_{3/2}(z0 exp(-y))/nu, x_bs=alpha exp(-gamma*y),
    (xi^2*y')'=xi^2*(x_n+x_bs). Implies u^2=xi*y'.
    Alpha=0 produces a purely noncondensed equilibrium.
    """
    z0=central_fugacity(nu)
    eps=1.e-5
    w0=[(1+alpha)*eps**2/6., (1+alpha)*eps**3/3.]
    def dydxi(xi, w):
        y,q=w
        z=z0*np.exp(-max(y,0.))
        density=float(mp.polylog(1.5,z))/nu
        bs=alpha*np.exp(-gamma*y)
        return [q/xi**2, xi**2*(density+bs)]
    xi=np.linspace(eps,xi_max,320)
    sol=solve_ivp(dydxi,[eps,xi_max],w0,t_eval=xi,
                  rtol=1e-8,atol=1e-11)
    if not sol.success:raise RuntimeError(sol.message)
    y,q=sol.y
    xn=np.array([float(mp.polylog(1.5,z0*np.exp(-max(t,0.))))/nu for t in y])
    xb=alpha*np.exp(-gamma*y)
    assert 0 <= xn.min() and xn.max() <= 1.0000001
    assert np.max(np.diff(xn)) < 1e-7
    return {"xi":xi.tolist(),"x_n":xn.tolist(),"x_bs":xb.tolist(),
            "u":np.sqrt(q/xi).tolist(),"z0":float(z0)}

def mass_eV(nu,R_kpc,V_km_s):
    m_kg=(np.sqrt(np.pi/2)*hbar**3 /
          (G*nu*(R_kpc*1000*parsec)**2*(V_km_s*1000)))**.25
    return m_kg/(elementary_charge/c**2)

def main():
    root=Path(__file__).resolve().parents[1]
    (root/"figures").mkdir(parents=True,exist_ok=True)
    (root/"code").mkdir(parents=True,exist_ok=True)
    nu_values=[0.03,0.5,1.5,2.2,NU_CRIT]
    prof={str(v):equilibrium(v) for v in nu_values}
    crit=equilibrium(NU_CRIT,xi_max=0.01)
    critical_slope=2*np.sqrt(np.pi)/(np.sqrt(6)*NU_CRIT)
    measured_slope=(1.-crit["x_n"][0])/crit["xi"][0]
    m=mass_eV(2.61,0.6,44.1)
    assert abs(measured_slope-critical_slope) < 2e-4
    assert abs(m-86.2245655) < 0.01
    audit={"nu_critical":NU_CRIT,"critical_slope_analytic":critical_slope,
           "critical_slope_numeric_at_xi_1e-5":measured_slope,
           "source_F583_4_mass_eV":m,
           "status":"equilibrium ODE and one source formula verified; no galaxy data fit"}
    (root/"code"/"audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for variable,filename,ylabel in [
        ("x_n","bose_density.png","rho_n/rho_n(0)"),
        ("u","bose_circular_velocity.png","V/V0")]:
        fig,ax=plt.subplots(figsize=(7,4.5))
        for nu,p in prof.items():
            ax.plot(p["xi"],p[variable],label="nu="+str(round(float(nu),3)))
        ax.set(xlabel="xi=r/R",ylabel=ylabel,xlim=(0,12))
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(root/"figures"/filename,dpi=140)
        plt.close(fig)
    print(json.dumps(audit,indent=2))

if __name__=="__main__":
    main()
