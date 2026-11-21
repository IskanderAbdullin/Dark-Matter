#!/usr/bin/env python3
"""Two-component Poisson–Bose synthetic comparison, not a galaxy fit."""
import json
from pathlib import Path
import numpy as np
from reproduce_profiles import equilibrium
NU_CRIT=2.612375348685488
alpha=gamma=12.8
root=Path(__file__).resolve().parents[1]
radii=[0.1,0.5,1.,3.,6.,10.]
record={}
for nu in [0.19,NU_CRIT]:
    p=equilibrium(nu,xi_max=12)
    t=equilibrium(nu,alpha=alpha,gamma=gamma,xi_max=12)
    record[str(nu)]=[]
    for xi in radii:
        pure=float(np.interp(xi,p["xi"],p["u"]))
        mixed=float(np.interp(xi,t["xi"],t["u"]))
        record[str(nu)].append({"xi":xi,"u_pure":pure,"u_two_component":mixed,
                                 "ratio":mixed/pure})
a=2*np.sqrt(np.pi)/(np.sqrt(6)*NU_CRIT)
near=equilibrium(NU_CRIT,alpha=alpha,gamma=gamma,xi_max=.01)
record["critical_slope_analytic"]=a*np.sqrt(1+alpha)
record["critical_slope_numeric_at_xi_1e-5"]=(1-near["x_n"][0])/near["xi"][0]
(root/"code"/"two_component_results.json").write_text(json.dumps(record,indent=2))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(8,5))
for nu in [0.19,NU_CRIT]:
    p=equilibrium(nu,xi_max=12);t=equilibrium(nu,alpha=alpha,gamma=gamma,xi_max=12)
    ax.plot(p["xi"],p["u"],label=f"pure nu={nu:.3g}")
    ax.plot(t["xi"],t["u"],"--",label=f"two-component nu={nu:.3g}")
ax.set(xlabel="xi=r/R",ylabel="V/V0",title="Synthetic equilibrium; not observed rotation curves",xlim=(0,12))
ax.legend();fig.tight_layout()
(root/"figures").mkdir(parents=True,exist_ok=True)
fig.savefig(root/"figures"/"two_component_comparison.png",dpi=160)
plt.close(fig)
print(json.dumps(record,indent=2))
