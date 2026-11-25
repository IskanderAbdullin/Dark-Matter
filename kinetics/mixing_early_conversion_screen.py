"""Conditional early heavy-state seeding/catalysis and late threshold tests.

The light-state PSD is held fixed; no claim of full coupled quantum
evolution or pointwise bound including arbitrary incoming heavy PSD.
Conservative AA gain bound ignores thresholds and uses E_H>=m and
angular interval <=2 p_H. The extra factor 4 is a normalization margin.
Catalytic LH->HH count bounds number growth if f_H<=1.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from quartic_consistency import shell_nodes
from late_activation_overlap import H
from soft_mixing_activation import frac_B
from direct_activation_consistency import T0,CGS_RHO_TO_EV4
from resonant_mediator_inverse import HBAR,M
from velocity_only_exposure import RHO0

def analyze(row,case):
    p=case['params'];hist=case['result']['history']
    r=np.array([h['z']+1 for h in hist]);x=np.array([h['x'] for h in hist])
    # Remove any duplicated z=0 final sample.
    _,ix=np.unique(r,return_index=True);ix=np.sort(ix);r,x=r[ix],x[ix]
    u=frac_B(x,p['delta'],p['Delta']);t=-np.log(r)
    theta=row['final']['theta'];Tref=T0*(3.909/row['g_bath'])**(1/3)/theta
    e=np.asarray(row['frozen_cell_spectrum']['edges'])
    f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.);fmax=float(f.max())
    P=Tref*r[:,None]*e[None,:]
    E=np.hypot(P,M)
    dE=(P[:,1:]**2-P[:,:-1]**2)/(E[:,1:]+E[:,:-1])
    I=dE@f
    q,w=shell_nodes(e,f,12);w/=w.sum()
    pphys=Tref*r[:,None]*q[None,:];Ep=np.hypot(pphys,M)
    n=RHO0*CGS_RHO_TO_EV4/M*r**3
    variance=n*np.sum(w[None,:]/Ep,axis=1)
    hubble=HBAR*H(r-1)
    lam=p['lambda_B'];ampHH2=lam**2*u*u*(1-u)**2
    ampLH2=lam**2*u**3*(1-u)
    ampCat2=lam**2*u*(1-u)**3
    # Conservative per-mode gain on f_H=0, including identical/mixed
    # final state convention margin. If f_H<=1, gain multipliers <=4.
    gain= (ampHH2+2*ampLH2*(1+fmax))*I**2/(32*np.pi**3*M*hubble)
    # For one H with E_H>=m, cross-section times relative speed is bounded
    # by amplitude^2/(32 pi E_H E_A) for the contact process. Add margin.
    cat=4*ampCat2*variance/(32*np.pi*M*hubble)
    # AA->HH and AA->LH number gain relative to light number, conservative
    # f_H<=1 estimate; extra factor 4 allows both produced heavy legs.
    source=4*(ampHH2+2*ampLH2*(1+fmax))*variance**2/(32*np.pi*n*hubble)
    early=r>=1001
    integral=lambda y:float(trapezoid(y[early],t[early]))
    es=integral(source);ec=integral(cat)
    # Only under fixed reservoir and bounded heavy occupation; cannot
    # infer f_H<=1 merely from this global number estimate.
    logbound=np.log(max(es,1e-300))+ec
    kinetic=np.sum(w[None,:]*pphys**2/(Ep+M),axis=1)
    min_gap=p['Delta']/(2*M)
    zrows=[]
    for z in [1.,3.,5.,10.,100.]:
        j=np.argmin(abs(np.log(r/(1+z))))
        zrows.append(dict(z=float(r[j]-1),kinetic_energy_per_particle_eV=float(kinetic[j]),
            instantaneous_heavy_fraction_energy_cap=float(kinetic[j]/min_gap)))
    # Adiabaticity diagnostic in late regime, relative to vacuum level
    # splitting; excludes possible early finite-temperature level crossing.
    dxdt=np.gradient(x,t)
    angdot=abs(p['delta']*p['Delta']*dxdt/(p['Delta']**2+4*p['delta']**2*x*x))*hubble
    gap=np.sqrt(p['Delta']**2+4*p['delta']**2*x*x)/(2*M)
    return dict(seed=case['result']['seed'],muH=p['mu_over_H_on'],fA_max=fmax,
        integrated_AA_heavy_mode_seed_upper_frozen=integral(gain),
        integrated_LH_to_HH_per_heavy_count_upper_if_fH_le_1=ec,
        integrated_AA_heavy_number_source_upper_if_fH_le_1=es,
        log_heavy_fraction_upper_conditional=min(logbound,1e300),
        early_peak_gain_per_H=float(max(gain[early])),
        early_peak_catalytic_per_H=float(max(cat[early])),
        late_vacuum_adiabaticity_max=float(max((angdot/gap)[r<=1001])),
        late_energy_caps=zrows,min_gap_eV=min_gap,
        scope='AA source on frozen A PSD; catalytic estimate conditional on fH<=1, no bootstrap for full momentum distribution. Instantaneous kinetic-energy cap at fixed gap, excludes field work. No full finite-temperature flavor kinetic solution.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--cases',required=True);a.add_argument('--output',required=True);x=a.parse_args()
    row=json.loads(Path(x.input).read_text())[0]
    out=[analyze(row,c) for c in json.loads(Path(x.cases).read_text())]
    Path(x.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in out:
        print({k:v for k,v in r.items() if k not in ['scope','late_energy_caps']},flush=True)
    print(out[-1]['late_energy_caps'])
