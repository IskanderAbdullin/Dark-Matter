"""Detuning-independent early flavor-transfer bound and postportal adiabatic test.

For a positive two-state density matrix |rho_AB|<=sqrt(rho_A rho_B).
With purely diagonal damping/source and no direct B source, B/N_injected
<= (integral |H_AB| dt)^2, independent of diagonal detuning and dephasing.
This includes arbitrary thermal crossings within the covered interval.
It is NOT a full collision operator or proof that quasiparticles exist.
After a specified portal-off epoch, assume vacuum mass matrix and use
the instantaneous-eigenbasis off-diagonal rotation. Integration by parts
bounds the late transition amplitude for a smooth given history.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad,solve_ivp,trapezoid
from soft_mixing_activation import scales,variance_table,bath_integrals
from late_activation_overlap import H,H0,Om,Or
from direct_activation_consistency import CGS_RHO_TO_EV4
from resonant_mediator_inverse import HBAR,M
from velocity_only_exposure import RHO0

def vacuum_energy(p,x,d,D):
    split=np.hypot(D,2*d*x)
    Eplus=np.sqrt(p*p+M*M+(D+split)/2)
    Eminus=np.sqrt(p*p+M*M+(D-split)/2)
    return split/(Eplus+Eminus)

def early_bound(row,p,seed=1e-8,rstart=1e16,rcut=1e14):
    # Early coarse expansion bound: H >= half of radiation-only LCDM
    # evaluated at same redshift; g*/gs corrections must respect this.
    V,Tref=variance_table(row);delta=p['delta'];D=p['Delta'];mu=p['mu_eV'];X=p['X']
    hc=.5*H0*np.sqrt(Or)*HBAR
    n0=RHO0*CGS_RHO_TO_EV4/M
    # If the actual A number never exceeds twice final comoving abundance,
    # <A^2> <= 2*n0*r^3/m. Disregard damping, keep absolute restoring force.
    # d(ln|chi|) bound via integral of |Q| with friction >=1 in radiation.
    b=mu**2/hc**2*(1/(4*rstart**4)+1/(2*(1+p['z_on'])**2*rstart**2)
        +2*delta**2*n0/(D*X*M*rstart))
    pre_factor=np.exp(b)
    seed_max=seed*pre_factor
    # Derivative portal-era exact positive energy square root:
    # |H_AB|=|delta*x|/(Eplus+Eminus) <= |delta*x|/sqrt(p^2+m^2)
    # if B diagonal mass >=m and positive total squared-energy matrix.
    # Use seed margin 1.01, verified by the homogeneous roll below.
    maxchi=seed_max*1.01
    a=np.array([0.,row['final']['peak_bin_q'],.1,1.,max(row['frozen_cell_spectrum']['edges'])])
    # Integrate r=rcut/u, u in [0,1] over the arbitrary-detuning era.
    vals=[]
    for q in a:
        I=abs(delta)*maxchi/(hc*rcut**2)*quad(lambda t:t*t/np.hypot(M*t,q*Tref*rcut),0.,1.,epsabs=1e-18,epsrel=1e-10)[0] if q>0 else abs(delta)*maxchi/(2*hc*rcut**2*M)
        vals.append(dict(q_comoving=float(q),amplitude_bound=float(I),probability_bound=float(I*I)))
    return dict(rcut=rcut,rstart=rstart,seed=seed,pre_roll_growth_exponent_bound=float(b),assumed_seed_envelope=maxchi,
       x_parent_at_rcut=86e12/(Tref*rcut),mode_bounds=vals,
       scope='Positive quasiparticle Hamiltonian, B diagonal mass >=m, H>=0.5 radiation extrapolation, no direct B injection, initial A trace <= total injection, no off-diagonal bath term. Early chi bound assumes number<=2 final and friction>=1.')

def solve_roll(row,p,scheme,seed=1e-8,rtol=1e-9):
    V,Tr=variance_table(row);eta=bath_integrals()
    D=p['Delta'];d=p['delta'];X=p['X'];mu=p['mu_eV'];z=p['z_on']
    eta0=p['m_eta_at_vev_eV'];td0=p['Td_on_eV']/(1+z)
    cw=p['Neta']*eta0**4/(64*np.pi**2)
    def rhs(t,y):
        r=np.exp(-t);x,dx=y;hh=H(r-1)*HBAR
        dh=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        therm=12*(r/(1+z))**2*eta(eta0*abs(x)/(td0*r))
        med=d*d*V(r)/(X*np.hypot(D,2*d*x))
        logx=np.log(max(x*x,1e-300))
        loop=4*cw/X*((x*x*logx) if scheme=='origin-vev' else (1+x*x*(logx-1)))
        return [dx,-(3+dh)*dx-(mu/hh)**2*x*(x*x-1+therm-med+loop)]
    sol=solve_ivp(rhs,(-np.log(1e16),0.),[seed,0.],rtol=rtol,atol=1e-13,
       method='DOP853',max_step=.03,dense_output=True)
    assert sol.success
    return sol,Tr

def post_bound(row,p,scheme,seed=1e-8,rcut=1e14):
    early=early_bound(row,p,seed,rcut=rcut)
    sol,Tr=solve_roll(row,p,scheme,seed)
    D=p['Delta'];d=p['delta'];qmax=max(row['frozen_cell_spectrum']['edges'])
    ti=-np.log(rcut);tm=-np.log(1001)
    # Until z=1000 use absolute variation of the vacuum angle; no
    # adiabatic assumption is needed for this interval.
    def vals(t):
        y=sol.sol(t);x=y[0];dx=y[1]
        ang=.5*np.arctan2(2*d*x,D)
        g=d*D*dx/(D*D+4*d*d*x*x)
        return ang,g,x
    t=np.linspace(-np.log(1e16),ti,200)
    mx=float(max(abs(sol.sol(t)[0])))
    assert mx<=early['assumed_seed_envelope']
    t=np.linspace(ti,tm,10000);ang,g,x=vals(t)
    TV=float(trapezoid(abs(g),t))
    # Initial heavy eigenmode amplitude: project initially A flavor
    # and add the early arbitrary-detuning bound.
    entry=abs(np.sin(ang[0]))+np.sqrt(early['mode_bounds'][0]['probability_bound'])+TV
    out=[]
    for q in [0.,row['final']['peak_bin_q'],qmax]:
        res=[]
        for nn in [4000,16000,64000]:
            t=np.linspace(tm,0.,nn);a,g,x=vals(t);r=np.exp(-t)
            W=vacuum_energy(Tr*q*r,x,d,D)/(HBAR*H(r-1))
            s=g/W
            # Integration by parts with |a_L|<=1, |a'_L|<=|g|.
            # Discrete total variation and trapezoid checked for convergence;
            # this is a numerical estimate of the analytic bound, not
            # an interval-arithmetic proof.
            ad=abs(s[0])+abs(s[-1])+np.sum(abs(np.diff(s)))+trapezoid(abs(s*g),t)
            res.append(dict(samples=nn,adiabatic_amplitude_bound_estimate=float(ad),
                max_abs_g_over_W=float(max(abs(s))),minimum_W=float(min(W))))
        out.append(dict(q=q,refinements=res,combined_heavy_eigenstate_probability_bound_estimate=float((entry+res[-1]['adiabatic_amplitude_bound_estimate'])**2)))
    return dict(early=early,scheme=scheme,max_pre_rcut_field=mx,vacuum_angle_at_rcut=float(ang[0]),
        pre_late_angle_total_variation=TV,pre_late_heavy_eigenstate_amplitude_bound=entry,
        late=out,today_field=float(sol.y[0,-1]),
        scope='Postcut vacuum diagonal mass difference assumed, no scattering/coherence injection. Numerical integral bound for specified homogeneous chi solution; thermal/density corrections after cutoff and direct B collision source not covered.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--cases',required=True);a.add_argument('--output',required=True);v=a.parse_args()
    row=json.loads(Path(v.input).read_text())[0]
    cases=json.loads(Path(v.cases).read_text());out=[]
    for case in cases:
        p=case['params'];r=post_bound(row,p,case['result']['loop_scheme'])
        out.append(dict(params=p,result=r))
        Path(v.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print(p['sigma'],r['early']['mode_bounds'],r['pre_late_angle_total_variation'],r['pre_late_heavy_eigenstate_amplitude_bound'],r['late'],flush=True)
