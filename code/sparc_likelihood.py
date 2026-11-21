#!/usr/bin/env python3
"""SPARC mass-model parser and a Bose halo forward likelihood (NOT an executed fit).

Official source https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt
Lelli, McGaugh & Schombert (2016). Vobs random errors exclude inclination
systematics, which require treatment in publishable galaxy population inference.
"""
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import argparse
import numpy as np
from scipy.constants import G,hbar,c,e as e_charge,parsec
from reproduce_profiles import equilibrium

NU_CRIT=2.612375348685488

@dataclass(frozen=True)
class Datum:
    name:str
    D:float
    R:float
    Vobs:float
    sigma:float
    Vgas:float
    Vdisk:float
    Vbul:float
    SBdisk:float
    SBbul:float

def read_sparc(path):
    obs=defaultdict(list)
    for line_number,line in enumerate(Path(path).read_text(encoding="utf8").splitlines(),1):
        parts=line.split()
        if len(parts)!=10:continue
        try: xs=[float(v) for v in parts[1:]]
        except ValueError:continue
        if xs[1]<=0 or xs[3]<=0:raise ValueError(f"Invalid R/sigma at line {line_number}")
        d=Datum(parts[0],*xs)
        obs[d.name].append(d)
    if not obs:raise ValueError("No data rows found")
    for g,rows in obs.items():
        rs=[x.R for x in rows]
        if len(set(rs))!=len(rs):raise ValueError(f"Repeated radius: {g}")
    return dict(obs)

def nu_from_mass(m_eV,R_kpc,V0_kms):
    mkg=m_eV*e_charge/c**2
    return np.sqrt(np.pi/2)*hbar**3/(G*mkg**4*(R_kpc*1000*parsec)**2*(V0_kms*1000))

def predict_model(galaxy,m_eV,R_kpc,V0_kms,disk_ML=0.5,bulge_ML=0.7):
    if min(m_eV,R_kpc,V0_kms)<=0:raise ValueError("positive mass, R and V0 required")
    if min(disk_ML,bulge_ML)<0:raise ValueError("nonnegative mass-to-light required")
    nu=nu_from_mass(m_eV,R_kpc,V0_kms)
    if not 0<nu<=NU_CRIT:raise ValueError("unphysical central Bose degeneracy")
    xi=np.array([p.R/R_kpc for p in galaxy])
    profile=equilibrium(float(nu),xi_max=max(20.,1.1*max(xi)))
    u=np.interp(xi,profile["xi"],profile["u"])
    # SPARC Vgas can be negative; gas's contribution is signed Vgas*|Vgas|.
    gas=np.array([p.Vgas*abs(p.Vgas) for p in galaxy])
    vd=np.array([p.Vdisk for p in galaxy]);vb=np.array([p.Vbul for p in galaxy])
    V2=V0_kms**2*u**2+gas+disk_ML*vd**2+bulge_ML*vb**2
    if np.any(V2<0):raise ValueError("negative predicted V^2")
    return np.sqrt(V2)

def diagonal_loglike(galaxy,m_eV,R_kpc,V0_kms,disk_ML=0.5,bulge_ML=0.7):
    v=predict_model(galaxy,m_eV,R_kpc,V0_kms,disk_ML,bulge_ML)
    measured=np.array([p.Vobs for p in galaxy]);sigma=np.array([p.sigma for p in galaxy])
    return float(-.5*np.sum(((measured-v)/sigma)**2+np.log(2*np.pi*sigma**2)))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("official_sparc_table",type=Path)
    p.add_argument("--galaxy")
    p.add_argument("--mass-eV",type=float,default=86.0)
    p.add_argument("--R-kpc",type=float,default=0.6)
    p.add_argument("--V0-kms",type=float,default=44.1)
    args=p.parse_args()
    data=read_sparc(args.official_sparc_table)
    print("Galaxies",len(data),"data points",sum(map(len,data.values())))
    if args.galaxy:
        if args.galaxy not in data:raise ValueError("Unknown galaxy")
        print("Unfitted fixed-parameter diagonal log-likelihood:",
            diagonal_loglike(data[args.galaxy],args.mass_eV,args.R_kpc,args.V0_kms))
        print("This is NOT a fitted result or publication-ready likelihood.")
if __name__=="__main__":main()
