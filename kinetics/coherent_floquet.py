"""Conditional coherent real-field self-resonance for V=m^2 phi^2/2+lambda phi^4/24.

Homogeneous flat-space oscillator at one cosmological snapshot. Jacobi cn
is exact, fluctuation equation is n=2 Lame. Linear instability is not
proof that the incoherent production spectrum actually forms this field,
nor does it determine the final nonlinear spectrum.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import ellipj,ellipk
from scipy.integrate import solve_ivp,quad
from scipy.optimize import minimize_scalar
from quartic_consistency import MPL

def band(s):
    """Elliptic parameter s=k_elliptic^2 in [0,1/2]; return kappa^2 interval."""
    lo=3-3*s
    width=3*s*s/(2*np.sqrt(1-s+s*s)+2-s)
    return lo,lo+width

def floquet(s,k2,rtol=2e-11):
    alpha=1-2*s;beta=2*s
    period=2*ellipk(s) # period of cn^2, half the field period
    def rhs(t,z):
        y=ellipj(t,s)[1];omega2=k2+alpha+3*beta*y*y
        return [z[1],-omega2*z[0],z[3],-omega2*z[2]]
    sol=solve_ivp(rhs,[0,period],[1,0,0,1],rtol=rtol,atol=rtol*.1,method='DOP853',max_step=.12)
    a,b,c,d=sol.y[:,-1]
    tr=(a+d)/2;det=a*d-b*c
    return dict(mu=float(np.arccosh(max(1.,abs(tr)))/period),trace_half=float(tr),det=float(det))

def peak(s,rtol=2e-11):
    lo,hi=band(s)
    if hi-lo<1e-9:return dict(k2=(lo+hi)/2,mu=0.,det=1.,unresolved_narrow_band=True)
    res=minimize_scalar(lambda p:-floquet(s,p,rtol)['mu'],bounds=(lo,hi),method='bounded',
                       options={'xatol':max((hi-lo)*1e-5,1e-11)})
    return dict(k2=float(res.x),**floquet(s,res.x,rtol),band_low=lo,band_high=hi)

def oscillator(U,m,lam):
    # Solve U = B Phi^2 + A Phi^4 with rationalized quadratic.
    A=lam/24;B=m*m/2
    Phi2=2*U/(B+np.sqrt(B*B+4*A*U))
    w=np.sqrt(m*m+lam*Phi2/6)
    s=lam*Phi2/(12*w*w)
    # J=(1/(2pi))*oint pi_phi dphi; action density in Tref^3.
    J=2/np.pi*np.sqrt(Phi2)*quad(lambda u:np.sqrt(max(0.,2*(U-B*Phi2*u*u-A*Phi2**2*u**4))),0,1,epsabs=1e-16,epsrel=1e-10)[0]
    return dict(U=U,Phi2=float(Phi2),frequency_scale=float(w),elliptic_parameter=float(s),action=float(J))

def screen(row,limits):
    out=[]
    for case in limits:
        lam=case['lambda_quartic'];states=[]
        for z in case['states']:
            if not z['equilibrium']['saturated']:continue
            x=z['x'];h=min(row['history'],key=lambda h:abs(h['x']-x))
            N=z['initial']['N'];U=z['initial']['U'];m=1e-12*x
            TH=MPL*x/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV'])*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
            for fraction in (.01,.1,1.):
                b=oscillator(U*fraction,m,lam)
                p=peak(b['elliptic_parameter']);lo,hi=band(b['elliptic_parameter'])
                omega=b['frequency_scale'];muH=p['mu']*omega*TH
                # Vacuum-level seeds in the linear resonance band; no actual gas matching.
                qa,qb=np.sqrt(lo)*omega,np.sqrt(hi)*omega
                mass_inst2=m*m+lam*b['Phi2']/2
                seed=quad(lambda q:q*q*np.sqrt(q*q+mass_inst2)/(4*np.pi**2),qa,qb,epsabs=1e-40)[0]
                # If each resonant mode grew with peak rate, this would be an optimistic
                # time to reach comparable energy, not a nonlinear calculation.
                efold=max(0.,.5*np.log(b['U']/max(seed,1e-300)))
                states.append(dict(x=x,energy_fraction=fraction,**b,**p,action_over_original_N=b['action']/N,
                    q_peak=np.sqrt(p['k2'])*omega,mu_over_H=muH,amplitude_efolds_from_vacuum=efold,
                    optimistic_H_time=efold/max(muH,1e-300),instantaneous_band_vacuum_energy_fraction=seed/b['U']))
        out.append(dict(sigma=case['sigma'],states=states))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--limits',required=True);p.add_argument('--output',required=True)
    a=p.parse_args()
    out=screen(json.loads(Path(a.input).read_text())[0],json.loads(Path(a.limits).read_text()))
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('quartic benchmark',peak(.5),flush=True)
    for c in out:
        print('SIGMA',c['sigma'],[(z['x'],z['energy_fraction'],z['mu_over_H'],z['optimistic_H_time']) for z in c['states'] if z['energy_fraction']==1.],flush=True)
