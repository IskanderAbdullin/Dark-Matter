"""Audit unrenormalized native CLASS output against exact archive moments."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points,snapshot_functions
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp


def native(path):
    with open(path) as fd:
        for line in fd:
            if "rho_ncdm[0]" in line:
                break
        else:
            raise ValueError("no ncdm columns")
    lookup={}
    for j,label in re.findall(r"(\d+):(.*?)(?=\s+\d+:|$)",line):
        lookup[label.strip()]=int(j)-1
    keys=["z","H [1/Mpc]","(.)rho_ncdm[0]","(.)p_ncdm[0]","(.)rho_ncdm[1]","(.)p_ncdm[1]"]
    a=np.loadtxt(path,usecols=[lookup[k] for k in keys])[::-1]
    assert np.all(np.diff(a[:,0])>0.)
    return a


def audit(table,shape,paths):
    bg=Background(table)
    s=json.loads(Path(shape).read_text())["snapshots"][-1]
    ai=float(bg.scale(s["T_visible_GeV"]))
    fp,_,fc,_=snapshot_functions(s)
    p=np.array(s["lab_momentum_edges_GeV"])
    s0=2*np.pi**2/45*bg.gs0*bg.T0**3
    factor=conversion(table)["Omega_h2_per_GeV_times_yield"]
    z=np.array([0.,1.,1000.,1.e5,1.e7,1.e8,1.e9,1./ai-1])
    xi=np.log1p(z)
    exact=[]
    native_tables=[native(path) for path in paths]
    rows=[]
    for species,g,mass,f,c in (("stable",4,m,fc,2),("scalar",1,mp,fp,4)):
        q,_,r=points(p,mass,f,48)
        qp=q*ai
        l=1+z[:,None]
        E=np.hypot(mass,l*qp)
        density=(E@r)*g*(ai*(1+z))**3
        pressure=(((l*qp)**2/(3*E))@r)*g*(ai*(1+z))**3
        omega=factor*density[0]/s0
        checks=[]
        for path,a in zip(paths,native_tables):
            xa=np.log1p(a[:,0])
            assert xa[0]<=xi[0] and xa[-1]>=xi[-1]
            rho=np.exp(PchipInterpolator(xa,np.log(a[:,c]))(xi))
            pr=np.exp(PchipInterpolator(xa,np.log(a[:,c+1]))(xi))
            native_omega=rho[0]/a[0,1]**2*bg.H0**2/1.e4
            checks.append(dict(class_background=Path(path).name,
                native_omega_h2=native_omega,
                relative_omega_error=native_omega/omega-1.,
                max_residual_shape_density_error=float(np.max(abs((rho/rho[0])/(density/density[0])-1.))),
                max_relative_pressure_to_density_error=float(np.max(abs((pr/rho)/(pressure/density)-1.))),
                rho_ratio_at_z=(rho/rho[0]).tolist(),
                pressure_over_energy=(pr/rho).tolist()))
        exact.append(dict(species=species,from_shell_energy_today_omega_h2=omega,
                       rho_ratio_at_z=(density/density[0]).tolist(),
                       pressure_over_energy=(pressure/density).tolist()))
        rows.append(dict(species=species,checks=checks))
    # Geometric validity of long-wave initial conditions: even at the
    # anchor later than radiation production, probe wavelengths are large.
    x,w=np.polynomial.legendre.leggauss(96)
    a=ai*(1+x)/2
    eta=float(bg.D0*ai/2*np.sum(w/np.sqrt(bg.weighted_radiation(a)+bg.om*a)))
    return dict(snapshot_absolute_a=ai,redshifts=z.tolist(),exact_shell_sources=exact,native_comparison=rows,
        approximate_conformal_horizon_Mpc_at_snapshot=eta,
        largest_k_h_per_Mpc_for_horizon_diagnostic=100.,
        k_times_approximate_conformal_horizon_at_100_h_per_Mpc=eta*bg.H0/100*100,
        scope="Background only; raw physical PSD with no input omega_ncdm rescaling. Earlier source perturbations, matter transfer, parameter likelihood, and exterior phase space are not audited here.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("shape");p.add_argument("class_backgrounds",nargs="+")
    a=p.parse_args()
    print(json.dumps(audit(a.table,a.shape,a.class_backgrounds),indent=2,allow_nan=False))
