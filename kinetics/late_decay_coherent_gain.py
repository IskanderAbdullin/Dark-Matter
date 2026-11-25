"""Undepleted-pump coherent small-signal toy for A->B+eta.
NR symmetric Doppler detunings, no expansion, B, medium or collisions.
The weak collective limit matches invariant occupation gain 2s.
This is a linear dispersion diagnostic, NOT a nonlinear field simulation.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from late_decay_bose_screen import mesh
from late_decay_line_gain import eta_gain
from quartic_consistency import shell_nodes,lambda_from_sigma
from direct_activation_consistency import T0
from resonant_mediator_inverse import M,HBAR,SEC_GYR

def dispersion(Omega,d,w):
    # Solve 1 = Omega^2 <atan(d/s)/(s*d)> for s.
    d=np.asarray(d);w=np.asarray(w)/sum(w)
    fun=lambda x: np.sum(w*np.where(d>0, np.arctan(d/(Omega*x))/(np.maximum(d,1e-300)/Omega*x),1/x**2))-1
    return Omega*brentq(fun,1e-24,1.0000001,xtol=1e-15,rtol=1e-13)

def run(row,sigma=.3,lifetime=10.):
    p0=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']*row['final']['peak_bin_q']
    e,_,f,meta=mesh(row,160,20.,0.)
    q,w=shell_nodes(e,f,12);N=w.sum();lam=float(lambda_from_sigma(sigma));gamma=1/(lifetime*SEC_GYR)
    out=[]
    for r in [1e4,1e5,1e6,2.5e6]:
        P=p0*r;n=N*P**3;v=P*q/M
        for ps in [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,.01]:
            g2=8*np.pi*M*M*gamma*HBAR/ps
            Om=np.sqrt(g2*n/(8*M*M*ps));dop=ps*v
            s=dispersion(Om,dop,w)
            kin=eta_gain(ps,r,ps,M,gamma,e,f,p0)
            nr= np.pi*Om**2*np.sum(w/dop)/N/HBAR
            mu=lam*n/(8*M*M)
            out.append(dict(r=r,recoil_eV=ps,pump_density_eV3=n,Omega_eV=Om,
                Omega_over_peak_doppler=Om/(ps*P/M),growth_occupation_per_s=2*s/HBAR,
                exact_invariant_kinetic_gain_per_s=float(kin),nr_weak_collective_gain_per_s=nr,
                coherent_to_NR_kinetic=2*s/HBAR/nr,meanfield_fullB_over_peak_doppler=mu/(ps*P/M),
                meanfield_fullB_over_recoil=mu/ps))
    return dict(scope='Linearized bosonic transition polarization with prescribed inverted gas, symmetric NR Doppler detuning and no phase-space backreaction, damping, medium shifts or redshift. Weak limit matches on-shell kinetic gain; strong regime flags failure of Markovian incoherent approximation.',cases=out,mesh=meta)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    v=run(json.loads(Path(a.input).read_text())[0]);Path(a.output).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    for x in v['cases']:
        if x['r']==1e5:print(x['recoil_eV'],x['Omega_over_peak_doppler'],x['coherent_to_NR_kinetic'],x['growth_occupation_per_s'])
