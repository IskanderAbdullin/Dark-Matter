"""Initial-condition and UV-scattering diagnostics, not a reheating model.

For a spectator B with quartic+quadratic potential, integrate the coherent
radiation-era onset and use the oscillator action after adiabaticity.
Fragmentation, inflaton couplings and isocurvature transfer are not included.
A derivative-portal contact freeze-in estimate uses massless MB bath
particles; it is not a full high-temperature reaction network.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.special import gamma
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
from derivative_portal_consistency import MP,MPL,GSTAR,M,R,K,m,RHO_OVER_S_DM,rate_mapping
from soft_mixing_activation import scales

def action(E,alpha,beta):
    if beta==0:A=np.sqrt(2*E/alpha)
    else:A=np.sqrt(4*E/(alpha+np.sqrt(alpha*alpha+4*beta*E)))
    def f(x):
        b=A*x
        return np.sqrt(max(0.,2*(E-.5*alpha*b*b-.25*beta*b**4)))
    return 2*A/np.pi*quad(f,0.,1.,epsabs=1e-20,epsrel=2e-9)[0]

def coherent(B,lam,xend=6000.,rtol=1e-8):
    w0=np.sqrt(m*m+lam*B*B/6);alpha=(m/w0)**2;beta=1-alpha
    t=np.linspace(xend*.9,xend,121)
    fun=lambda x,y:[y[1],-1.5*y[1]/x-alpha*y[0]-beta*y[0]**3]
    x0=1e-4
    sol=solve_ivp(fun,(x0,xend),[1-x0*x0/5,-2*x0/5],t_eval=t,rtol=rtol,atol=rtol*1e-5,max_step=10.)
    assert sol.success
    q=sol.y[0];E=.5*sol.y[1]**2+.5*alpha*q*q+.25*beta*q**4
    J=np.array([action(e,alpha,beta)*x**1.5 for e,x in zip(E,t)])
    Tscale=np.sqrt(w0*MP/(3.32*np.sqrt(GSTAR)));s=2*np.pi**2*GSTAR*Tscale**3/45
    frac=m*w0*B*B*float(J.mean())/s/RHO_OVER_S_DM
    return dict(B_initial_GeV=B,lambda_B=lam,omega0_GeV=w0,quadratic_fraction=alpha,
       fraction_DM_coherent_action=frac,adiabatic_action_relative_span=float((J.max()-J.min())/J.mean()),
       xend=xend,rtol=rtol)

def scattering():
    p=rate_mapping(M,R,K,.03,47.,GSTAR);F=2*p['f_over_C_GeV']
    # Contact from sterile field redefinition, one Weyl+scalar initial
    # component, spin sum, massless final state. Crossed channels and
    # quantum-statistical corrections may change the coefficient.
    sig=(p['y1']**2+p['y2']**2)/(32*np.pi*F*F)
    coef=4*sig*MP/(np.pi**4*(2*np.pi**2*GSTAR/45)*1.66*np.sqrt(GSTAR))
    Y=RHO_OVER_S_DM/m
    return dict(F_GeV=F,sigma_MB_GeV_minus2=sig,contact_Y_per_GeV=coef,
       required_Y=Y,rows=[dict(T_R_GeV=T,Y_MB=coef*T,fraction_DM_MB=coef*T/Y,
          fraction_with_illustrative_factor12=12*coef*T/Y) for T in [M*10,1e8,1e9,F/10,F]],
       scope='Massless MB one-contact benchmark with 4 EW/charge channels. Factor12 sensitivity, not a rigorous total upper bound; extrapolation to T_R=F is edge of EFT.')

def screen():
    p=scales(.1,z_on=5.,dneff=.1,muH=100.)
    lam=p['lambda_B'];ls=lam/6
    stochastic=(3/(2*np.pi**2*ls))**.25*np.sqrt(gamma(.75)/gamma(.25))
    rw=np.sqrt(60)/(2*np.pi)
    vals=[]
    # In quartic-dominated regime the dimensionless onset is universal.
    # Evaluate at representative amplitudes, check scaling and resolution.
    for B in [1e3,1e6,1e9,1e12]:
        vals.append(coherent(B,lam,xend=6000.))
    conv=coherent(1e9,lam,xend=24000.,rtol=2e-9)
    ref=vals[2];B1=ref['B_initial_GeV']*(.01/ref['fraction_DM_coherent_action'])**(2/3)
    Bverify=coherent(B1,lam,xend=6000.)
    v=p['vev_eV']*1e-9
    return dict(params=p,scattering=scattering(),coherent_spectator=vals,refinement=conv,
       B_at_1percent_scaling_GeV=B1,B_threshold_verification=Bverify,
       equilibrium_Brms_over_Hinf=float(stochastic),finite_60_efold_random_walk_Brms_over_Hinf=rw,
       Hinf_for_1percent_B_equilibrium_GeV=B1/stochastic,Hinf_for_1percent_B_random_walk_GeV=B1/rw,
       Hinf_for_chi_fluctuation_equal_seed_per_log_GeV=2*np.pi*v*1e-8,
       Hinf_for_chi_60efold_rw_equal_seed_GeV=v*1e-8/rw,
       instantaneous_reheat_H_min_GeV=1.66*np.sqrt(GSTAR)*(10*M)**2/MP,
       scope='Conditional radiation-era coherent action; quartic fragmentation and inflation/reheating model unspecified. Stochastic equilibrium requires long inflation; 60 e-fold free random walk assumes initially zero and negligible drift.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=screen()
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps(r,indent=2))
