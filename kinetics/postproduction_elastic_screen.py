"""Post-e+e- annihilation collisionless-PSD consistency and ideal equilibrium limits.

The frozen PSD is redshifted without alteration. Calculate the dilute elastic
collision count for a positive contact quartic using exact massive Moller
velocity and identical-particle cross section. If large, this is evidence
against assuming that collisions can be ignored, NOT the actual collisional
solution or a net thermalization rate.

Also compute an ideal number-conserving Bose-gas equilibrium limit at the
production endpoint and adiabatically evolve its saturated thermal component.
A real scalar has no exact conserved particle number; number-changing
processes and condensate dynamics are omitted deliberately.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from scipy.optimize import brentq
from quartic_consistency import shell_nodes,lambda_from_sigma,gap,variance

HBAR_GEV_S=6.582119569e-25
MPC_KM=3.0856775814913673e19
C_KMS=299792.458
T0_GEV=2.7255*8.617333262e-14
M=86e-9

def pair_factor(q,w,mu,angular_order=12):
    E=np.hypot(q,mu);sel=w>0
    q=q[sel];v=w[sel]/E[sel]
    if mu==0:return float(v.sum()**2)
    rapid=np.arcsinh(q/mu)
    d0=2*mu*mu*np.sinh((rapid[:,None]-rapid[None,:])/2)**2
    pp=q[:,None]*q[None,:]
    a,b=roots_legendre(angular_order)
    pair=0.
    for u,g in zip(a,b):
        d=d0+pp*(1-u)
        pair+=float(g*.5*(v@np.sqrt(d/(d+2*mu*mu))@v))
    return pair

def equilibrium_integrals(T,m):
    # Dimensionless p/T, chemical potential = rest mass (saturation).
    if T<=0:return (0.,0.,0.,0.)
    z=m/T
    def values(y):
        E=np.hypot(y,z)
        K=y*y/(E+z)
        if K>700:return np.zeros(4)
        f=1/np.expm1(K)
        # n, kinetic energy, pressure, entropy=(rho+P-m*n)/T
        return np.array([y*y*f,y*y*K*f,y**4*f/(3*E),
                         y*y*K*f+y**4*f/(3*E)])
    # NR thermal momenta scale as sqrt(z); choose finite upper cutoff.
    upper=max(70.,np.sqrt(140*z+4900.))
    from scipy.integrate import quad_vec
    u=quad_vec(values,0,upper,epsrel=2e-9,epsabs=1e-12)[0]/(2*np.pi**2)
    return float(T**3*u[0]),float(T**4*u[1]),float(T**4*u[2]),float(T**3*u[3])

def equilibrium_limit(row,q,w):
    mu=row['scalar_ratio']*row['xend']
    N=float(w.sum())
    kinetic=float(w@(q*q/(np.hypot(q,mu)+mu)))
    T=brentq(lambda t:equilibrium_integrals(t,mu)[1]-kinetic,1e-12,10.,xtol=1e-13)
    nth,uth,P,s=equilibrium_integrals(T,mu)
    if nth>=N:
        return dict(saturated=False,nth_over_N=nth/N,reason='No condensate implied by number and energy at endpoint')
    # Late NR saturated ideal limit: entropy per excited particle.
    zeta32=2.612375348685488;zeta52=1.341487257250917
    s_over_n_nr=2.5*zeta52/zeta32
    return dict(saturated=True,temperature_over_Tref=T,temperature_over_bath=T/row['final']['theta'],
        excited_number_fraction=nth/N,condensed_number_fraction=1-nth/N,
        entropy_per_total_particle=s/N,late_NR_excited_fraction=(s/N)/s_over_n_nr,
        kinetic_energy_over_Tref4=kinetic,pressure_over_Tref4=P,
        scope='Ideal saturation at fixed total number and energy, then adiabatic number-conserving elastic evolution; no Hartree/number-changing/condensate formation calculation.')

def screen(row):
    assert row['xend']>=60
    spec=row['frozen_cell_spectrum'];theta=row['final']['theta']
    q,w=shell_nodes(np.asarray(spec['edges']),np.asarray(spec['f_phi']),2)
    q4,w4=shell_nodes(np.asarray(spec['edges']),np.asarray(spec['f_phi']),4)
    N=float(w.sum());assert abs(N/row['final']['number']-1)<1e-8
    eq=equilibrium_limit(row,q4,w4)
    # Post annihilation (T_gamma<=30 keV), use standard radiation background.
    tr=(3.909/106.75)**(1/3);Tref0=T0_GEV*tr/theta
    amin=T0_GEV/(30e-6)
    aa=np.geomspace(amin,1.,72)
    h=.674;Ob=.02237/h**2;Od=row['final']['omega_phi_h2_m86']/h**2
    Or=2.4728e-5/h**2*(1+.22710731766*3.044);Ol=1-Ob-Od-Or
    H0=100*h/MPC_KM
    E0=float(w@np.hypot(M,q*Tref0)/N)
    unit_lam=float(lambda_from_sigma(1.))
    states=[]
    for a in aa:
        tref=Tref0/a;mu=M/tref
        Emean=float(w@np.hypot(M,q*tref)/N)
        H=H0*np.sqrt(Or/a**4+Ob/a**3+Od/a**3*Emean/E0+Ol)
        pf=pair_factor(q,w,mu,8)
        rate=unit_lam**2/(64*np.pi)*tref*pf/N/HBAR_GEV_S/H
        # Hartree screen for upper tested cross section at the same fixed PSD.
        u=gap(q4,w4,mu,float(lambda_from_sigma(1.25)))
        states.append(dict(a=float(a),z=float(1/a-1),Tgamma_eV=T0_GEV/a*1e9,
            H_per_s=float(H),horizon_k_per_Mpc=float(a*H*MPC_KM/C_KMS),
            rate_over_H_sigma1=float(rate),sigma_crit_rate_H_cm2_g=float(1/rate),
            max_test_gap_over_vacuum_mass=float(u/mu),
            mean_energy_over_mass=Emean/M))
    rates=np.array([s['rate_over_H_sigma1'] for s in states])
    ks=np.array([s['horizon_k_per_Mpc'] for s in states])
    horizons={}
    for k in [1.,10.,20.]:
        mask=aa<.01
        loga=float(np.interp(np.log(k),np.log(ks[mask][::-1]),np.log(aa[mask][::-1])))
        horizons[str(k)]=dict(z=float(np.exp(-loga)-1),
            sigma_crit_rate_H_cm2_g=float(np.exp(-np.interp(loga,np.log(aa),np.log(rates)))))
    cases=[]
    for s in [.001,.01,.1,1.25]:
        r=s*rates;idx=np.flatnonzero((r[:-1]>=1)&(r[1:]<1))
        zdec=None;kdec=None
        if len(idx):
            i=idx[-1];v=np.interp(0,np.log(r[i:i+2])[::-1],np.log(aa[i:i+2])[::-1])
            zdec=float(np.exp(-v)-1);kdec=float(np.exp(np.interp(v,np.log(aa),np.log(ks))))
        cases.append(dict(sigma_over_m_cm2_g=s,max_rate_over_H=float(r.max()),
            optical_depth_dilute=float(np.trapz(r,np.log(aa))),last_rate_H_crossing_z=zdec,
            horizon_k_at_crossing_per_Mpc=kdec))
    # Convergence at earliest sampled time (largest rate and most relativistic).
    mu=M/(Tref0/amin)
    d=abs(pair_factor(q4,w4,mu,16)/pair_factor(q,w,mu,8)-1)
    return dict(n=row['n'],kernel=row.get('kernel'),kappa=row['kappa'],omega=row['final']['omega_phi_h2_m86'],
        post_annihilation_start_Tgamma_keV=30.,states=states,cases=cases,
        horizon_entries=horizons,initial_pair_quadrature_drift=float(d),
        ideal_equilibrium_limit=eq,
        limitations='Vacuum-mass collisionless redshift, exact dilute elastic event count, no Bose enhancement/inverse cancellation. If rate is large, actual distribution will change. Standard radiation background after e+e- annihilation only. Ideal equilibrium limit assumes approximate number conservation, not a real-scalar formation proof.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();out=[]
    for f in a.input:
        row=json.loads(Path(f).read_text())[0];s=screen(row);s['input']=f;out.append(s)
        print(f,json.dumps(dict(cases=s['cases'],horizon_entries=s['horizon_entries'],ideal_equilibrium=s['ideal_equilibrium_limit'],
            initial_quadrature_drift=s['initial_pair_quadrature_drift'])),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
