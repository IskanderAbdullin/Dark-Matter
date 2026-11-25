"""Loop-scale splitting and linear pre-roll spatial stability diagnostics."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from soft_mixing_activation import scales,variance_table
from late_activation_overlap import H,Om,Or
from resonant_mediator_inverse import HBAR,HBARC,M
from direct_activation_consistency import MPL_REDUCED

def thresholds(p):
    u=p['u'];dmax=4*np.pi*np.sqrt(p['X_thermal_or_vacuum_cap'])*(1-2*u)/np.sqrt(u*(1-u))
    loop_split=p['lambda_B']*M*M/(32*np.pi**2)
    Lmax=np.sqrt(dmax*32*np.pi**2/p['lambda_B'])
    # For a rotating interaction basis protected by simultaneous O(2)
    # rotation, a radiatively generated aligned mass splitting makes the
    # inert eigenstate follow adiabatically under a Hubble-scale roll.
    gap=loop_split/(2*M);rate=p['mu_eV']
    return dict(max_split_eV2=dmax,one_loop_mass_scale_split_eV2=loop_split,
        max_unsuppressed_wilsonian_cutoff_eV=Lmax,cutoff_over_m=Lmax/M,
        aligned_basis_gap_eV=gap,homogeneous_rotation_rate_over_aligned_gap=rate/gap,
        maximum_relativistic_wall_thickness_m_for_nonadiabatic_rotation=HBARC/gap/100,
        scope='Illustrative Wilsonian quadratic threshold, no regulator-independent naturalness theorem. Exact symmetry rotates the induced splitting with interaction axis; adiabatic inert state then does not activate.')

def linear_modes(row,p,scheme='vev-curvature',seed=1e-8):
    V,_=variance_table(row);mu=p['mu_eV'];X=p['X'];delta=p['delta'];D=p['Delta'];z_on=p['z_on']
    acw=p['Neta']*p['m_eta_at_vev_eV']**4/(64*np.pi**2)
    cw0=4*acw/X if scheme=='vev-curvature' else 0.
    def eff(t):
        r=np.exp(-t);hh=HBAR*H(r-1)
        mass2=mu*mu*((r/(1+z_on))**2-1-delta*delta*V(r)/(D*X)+cw0)
        dh=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        return r,hh,dh,mass2
    def rhs(k):
        def f(t,y):
            r,hh,dh,ms=eff(t);kphys=k*r*HBARC/3.085677581e24
            return [y[1],-(3+dh)*y[1]-(ms+kphys*kphys)/hh**2*y[0]]
        return f
    t0=-np.log(1e8)
    def hit(t,y):return abs(y[0])-1e-3/seed
    hit.terminal=True
    zero=solve_ivp(rhs(0.),(t0,0.),[1.,0.],events=hit,rtol=1e-8,atol=1e-10,method='DOP853',max_step=.035,dense_output=True)
    end=zero.t[-1];zz=zero.y[0,-1];grid=np.linspace(t0,end,1800)
    band=[]
    for t in grid:
        r,hh,dh,ms=eff(t)
        band.append(np.sqrt(max(-ms,0.))/r*3.085677581e24/HBARC)
    out=[]
    for k in [0.,.0001,.001,.003,.01,.03,.1]:
        a=solve_ivp(rhs(k),(t0,end),[1.,0.],t_eval=grid,rtol=1e-8,atol=1e-10,method='DOP853',max_step=.035)
        assert a.success
        out.append(dict(k_Mpc_inv=k,final_transfer=float(a.y[0,-1]),relative_to_homogeneous_at_roll=float(abs(a.y[0,-1]/zz)),
                        maximum_absolute_transfer=float(max(abs(a.y[0])))))
    return dict(scheme=scheme,seed=seed,linear_stop_reached=bool(len(zero.t_events[0])),z_at_stop=float(np.exp(-end)-1),
       homogeneous_amplification=float(zz),maximum_pre_roll_tachyonic_k_Mpc_inv=float(max(band)),modes=out,
       scope='Linear spectator perturbations around chi=0, common initial amplitudes at z=1e8; metric, dark bath perturbations, nonlinear roll and domain wall dynamics omitted.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--output',required=True);x=a.parse_args()
    row=json.loads(Path(x.input).read_text())[0]
    pp=[scales(s,z_on=z,dneff=.1,muH=100.) for s,z in [(.1,5.),(.1,8.),(.3,8.),(1.25,8.),(1.25,10.)]]
    p=scales(.1,z_on=5.,dneff=.1,muH=100.,X=1.05*pp[0]['X_loop_floor'])
    result=dict(thresholds=[dict(params=q,diagnostics=thresholds(q)) for q in pp],
                linear=[linear_modes(row,p,s) for s in ['origin-vev','vev-curvature']])
    Path(x.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in result['thresholds']:print(r['params']['sigma'],r['params']['z_on'],r['diagnostics'])
    for v in result['linear']:print(v)
