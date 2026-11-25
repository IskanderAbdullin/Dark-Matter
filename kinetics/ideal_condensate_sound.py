"""Ideal saturated number-conserving Bose fluid: conditional sound horizon.

Assume instantaneous ideal equilibration of the production endpoint at fixed
scalar number and energy, then perfectly adiabatic local equilibrium of
condensate plus thermal excitations. A single tightly coupled perfect fluid
is an assumption: actual condensate/excitations can have distinct velocities.
Neglect mean-field pressure, number change and anisotropic stress. Calculate
only after e+e- annihilation, and do not identify sound horizon with a
measured matter transfer cutoff.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

def moments(z):
    # p/T = sqrt(1+m/T)*v resolves both relativistic and NR integrands.
    scale=np.sqrt(1+z)
    def fun(v):
        y=scale*v;E=np.hypot(y,z);K=y*y/(E+z)
        if K>700:return np.zeros(3)
        f=1/np.expm1(K)
        n=y*y*f
        return scale*np.array([n,n*K,n*y*y/(3*E)])
    n,u,p=quad_vec(fun,0,70,epsrel=5e-9,epsabs=1e-12)[0]/(2*np.pi**2)
    # n/T^3, u/T^4, P/T^4; s/T^3=u/T^4+P/T^4.
    return n,u,p
def run(row,scr,ngrid=360):
    eq=scr['ideal_equilibrium_limit'];assert eq['saturated']
    entropy_per_N=eq['entropy_per_total_particle']
    m=86.;Tgamma0=2.7255*8.617333262e-5
    tr=(3.909/106.75)**(1/3)
    theta=row['final']['theta'];N=row['final']['number']
    n0=N*(Tgamma0*tr/theta)**3
    s0=n0*entropy_per_N
    z=np.geomspace(1e-11,1e15,ngrid);T=m/z
    vals=np.array([moments(v) for v in z])
    n=T**3*vals[:,0];u=T**4*vals[:,1];p=T**4*vals[:,2]
    s=(u+p)/T
    # Standard thermodynamic consistency: dP/dT = s for mu=m.
    lp=PchipInterpolator(np.log(T[::-1]),np.log(p[::-1]))
    ls=PchipInterpolator(np.log(T[::-1]),np.log(s[::-1]))
    inv=PchipInterpolator(np.log(s[::-1]),np.log(T[::-1]))
    ln=PchipInterpolator(np.log(T[::-1]),np.log(n[::-1]))
    lu=PchipInterpolator(np.log(T[::-1]),np.log(u[::-1]))
    amin=Tgamma0/30000.
    a=np.geomspace(amin,1.,5000)
    sreq=s0/a**3
    lT=inv(np.log(sreq));Td=np.exp(lT)
    nt=np.exp(ln(lT));ut=np.exp(lu(lT));pt=np.exp(lp(lT))
    ntot=n0/a**3;rho=m*ntot+ut
    # d rho/dT=(m n/s+T) ds/dT along the isentropic equilibrium branch.
    cs2=Td/((m/entropy_per_N+Td)*ls.derivative()(lT))
    therm_error=np.max(np.abs(pt/Td*lp.derivative()(lT)/sreq-1))
    assert therm_error<.002,therm_error
    assert np.all(cs2>0) and np.all(cs2<1/3+1e-4)
    assert np.max(nt/ntot)<1
    h=.674;H0=100*h;Ob=.02237/h**2;Od=row['final']['omega_phi_h2_m86']/h**2
    Or=2.4728e-5/h**2*(1+.22710731766*3.044);Ol=1-Ob-Od-Or
    H=H0*np.sqrt(Or/a**4+Ob/a**3+Od*rho/rho[-1]+Ol)
    kh=a*H/299792.458
    rh=cumulative_trapezoid(299792.458*np.sqrt(cs2)/(a*a*H),a,initial=0.)
    eq_a=Or/(Ob+Od)
    sample={}
    for redshift in (1e8,1e7,1e6,2e5,1e5,1e4,3400,1000,0):
        aa=1/(1+redshift)
        if aa<amin:continue
        j=np.argmin(abs(np.log(a/aa)))
        sample[str(redshift)]=dict(z=float(1/a[j]-1),Tdark_eV=float(Td[j]),excited_fraction=float(nt[j]/ntot[j]),
            w=float(pt[j]/rho[j]),cs_over_c=float(np.sqrt(cs2[j])),sound_horizon_Mpc=float(rh[j]),
            instantaneous_sound_horizon_k_per_Mpc=float(kh[j]/np.sqrt(cs2[j])))
    r_eq=float(np.interp(np.log(eq_a),np.log(a),rh))
    return dict(ngrid=ngrid,entropy_per_particle=entropy_per_N,max_thermodynamic_error=float(therm_error),
        excited_fraction_today=float(nt[-1]/ntot[-1]),sound_horizon_to_equality_Mpc=r_eq,
        reciprocal_sound_horizon_to_equality_per_Mpc=1/r_eq,
        a_eq_CDM_reference=eq_a,sound_horizon_to_today_Mpc=float(rh[-1]),samples=sample,
        scope='Ideal saturated single perfectly coupled fluid; no condensate two-fluid modes, mean-field or number-changing; post-e+e- annihilation contribution only. Not a transfer prediction.')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--screen',required=True)
    ap.add_argument('--output',required=True);a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0];scr=json.loads(Path(a.screen).read_text())[0]
    results=[run(row,scr,n) for n in (240,480)]
    d=abs(results[0]['sound_horizon_to_equality_Mpc']/results[1]['sound_horizon_to_equality_Mpc']-1)
    assert d<.003,d
    Path(a.output).write_text(json.dumps(dict(results=results,relative_grid_drift=d),indent=2)+'\n')
    print(json.dumps(dict(eq_sound=results[-1]['sound_horizon_to_equality_Mpc'],
           today_sound=results[-1]['sound_horizon_to_today_Mpc'],drift=d,samples=results[-1]['samples'])),flush=True)
if __name__=='__main__':main()
