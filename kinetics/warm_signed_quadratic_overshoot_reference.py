"""Conditional positive signed quadratic Bose debit when a fermion cap exceeds 1/2."""
import argparse
import json

import numpy as np
from scipy.optimize import brentq
from warm_signed_fermion_time_comparison_reference import read, setup
from warm_fermion_depletion_kernel_reference import geom
from three_species_decay_onsager_reference import m,M,mp,y

def inverse(d,e,ns,nx,end=2.):
    s=np.geomspace(4*d*d,end,ns+1)
    x=np.linspace(0.,1.,nx+1)
    zl,zh=x[:-1][None,:]**2,x[1:][None,:]**2
    z=np.clip(.5,zl,zh);lo,hi=s[:-1,None],s[1:,None]
    al=e*e+lo/2-d*d-mp*mp; ah=e*e+hi/2-d*d-mp*mp
    pl,ph=np.maximum(lo/4-d*d,0),hi/4-d*d
    kl,kh=lo/4-mp*mp,hi/4-mp*mp
    Ll=(d+e)*al-2*d*kl*zh;Lh=(d+e)*ah-2*d*kh*zl
    den=geom(lo,d,e)+4*pl*kl*(1-zh)
    amp=8*ph*y**4*(hi*kh**2*z*(1-z)+np.maximum(Ll**2,Lh**2))/den**2
    beta=np.sqrt(1.-4*d*d/hi[:,0])
    v=(amp@(2*np.diff(x)))*d*d*beta/(np.pi*lo[:,0]**2)
    return max(float(max(v)),4*d*d/end*y**4/(np.pi*np.sqrt(geom(end,d,e))))*(1+1e-11)

def lower(model,cap,c):
    W,K,J1,J2,Z,eta,cut=model
    Y=np.ones((len(eta),cut.size))
    a=cap**(1/3)
    J1=J1*a;J2=J2*a*a
    c=c*max(0,2*cap-1)
    for s,h,t in W:
        r=np.cbrt(Y)
        Y=np.maximum(0.,Y-s*K*Y*Y-h*J2*r**4-t*J1*r*r-s*c-2e-14)
    return max(0.,float(np.max(Y-Z)))


def run(fermion,leptons,isolated,reference):
    F,L,I,R=map(read,(fermion,leptons,isolated,reference))
    cells=(320,620)
    K=[[inverse(d,e,n,2*n)/d**2 for d,e in ((m,M),(M,m))] for n in cells]
    assert all(x<=y for x,y in zip(K[-1],K[0],strict=True))
    rows=[]
    for a,b,c,r in zip(F["cases"],L["cases"],I["cases"],R["cases"],strict=True):
        assert a["q"]==b["q"]==c["q"]==r["q"]
        n,o=a["core_start_fermion_number_GeV3"],a["conservative_initial_fermion_Omega"]
        E=b["available_energy_upper_GeV4"]
        C=E*E*sum(K[-1])/n
        model=setup(n,E,L["coefficient"],c["bounds"][1]["clock"],2800)
        vals=[dict(fermion_per_spin_cap=f,lower_Omega_h2=lower(model,f,C)*o)
              for f in (.49,.5,.50005,.5001,.5002,.5005,.501)]
        assert abs(vals[1]["lower_Omega_h2"]-r["controls"][-1]["fermion_Omega_lower_at_half_per_spin"])<1e-12
        z=brentq(lambda f:lower(model,f,C)*o-F["target"],.5,.6,xtol=1e-12)
        rows.append(dict(q=a["q"],extra_rate_coefficient_before_excess_factor=C,
                         last_certifying_cap_approx=z,caps=vals))
    return dict(scope="Conditional homogeneous isotropic signed vacuum network. Bounding "
                      "a possible positive quadratic Bose term yields only a very narrow "
                      "extension past the half-per-spin condition. No achieved cap inferred.",
                scalar_inverse_energy_weighted_coefficients_by_s_grid=K,s_cells=cells,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("fermion","leptons","isolated","reference"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.fermion,a.leptons,a.isolated,a.reference),indent=2,allow_nan=False))
