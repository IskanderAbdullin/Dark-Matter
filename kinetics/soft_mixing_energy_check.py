"""Independent field+bath+light-matter energy accounting for soft activation.

Subtract free massless eta and constant-m light matter backgrounds.
Evolve work variable along with homogeneous roll in fixed LCDM H.
First order in tiny A mass shift; exact eta energy and pressure quadrature.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,zeta
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from soft_mixing_activation import scales,frac_B,variance_table,bath_integrals
from direct_activation_consistency import T0,CGS_RHO_TO_EV4,rho_tot
from quartic_consistency import shell_nodes
from late_activation_overlap import H,Om,Or
from resonant_mediator_inverse import HBAR,M,C
from velocity_only_exposure import RHO0

def bath_energy_table():
    u,w=roots_legendre(48);e=np.array([0.,3e-5,3e-4,.003,.03,.3,1.,4.,12.,50.])
    y=((e[:-1,None]+e[1:,None])/2+np.diff(e)[:,None]*u/2).ravel()
    w=(np.diff(e)[:,None]*w/2).ravel()
    wt=w*y*y/np.expm1(y)/(2*np.pi**2)
    j=np.geomspace(1e-6,1e4,1000)
    E=np.hypot(y[None,:],j[:,None])
    D=np.sum(wt[None,:]*j[:,None]**2/(E+y[None,:]),axis=1)
    V=np.sum(wt[None,:]/E,axis=1)
    d=PchipInterpolator(np.log(j),np.log(D));v=PchipInterpolator(np.log(j),np.log(V))
    def f(s):
        if s<1e-6:return s*s/24,1/12
        if s>1e4:return s*zeta(3,1)/np.pi**2-np.pi**2/30,zeta(3,1)/(np.pi**2*s)
        return np.exp(d(np.log(s))),np.exp(v(np.log(s)))
    return f

def run(row,p,scheme,rtol,seed=1e-8):
    V,Tref=variance_table(row)
    e=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    q,w=shell_nodes(e,f,12);w/=w.sum()
    r=np.geomspace(1.,1e16,1200);mom=r[:,None]*Tref*q[None,:];E=np.hypot(M,mom)
    W=RHO0*CGS_RHO_TO_EV4/M*r**3*np.sum(w[None,:]*mom*mom/E**3,axis=1)
    ww=PchipInterpolator(np.log(r),np.log(np.maximum(W,1e-300)))
    eta=bath_energy_table();X=p['X'];mu=p['mu_eV'];v=p['vev_eV'];delta=p['delta'];D=p['Delta']
    m_eta=p['m_eta_at_vev_eV'];td0=p['Td_on_eV']/(1+p['z_on']);Neta=p['Neta']
    acw=Neta*m_eta**4/(64*np.pi**2)
    def pot(x):
        x2=x*x;lx=np.log(max(x2,1e-300))
        if scheme=='origin-vev':
            val=acw*(x2*x2*(lx-.5)+.5);force=4*acw*x2*lx
        elif scheme=='vev-curvature':
            val=acw*(x2*x2*(lx-1.5)+2*x2-.5);force=4*acw*(1+x2*(lx-1))
        else:val=force=0.
        return X/4*(x2-1)**2+val,X*(x2-1)+force
    def local(t,x,dx):
        r=np.exp(-t);hh=H(r-1)*HBAR;Td=td0*r;s=m_eta*abs(x)/Td
        dr,ev=eta(s)
        er=Neta*Td**4*dr;ep=Neta*Td**4*(dr-s*s*ev)/3
        dm2=-2*delta*delta*x*x/(np.hypot(D,2*delta*x)+D)
        ar=.5*dm2*V(r);ap=-dm2*np.exp(ww(np.log(np.clip(r,1,1e16))))/6
        U,F=pot(x);K=.5*v*v*hh*hh*dx*dx
        force=F+Neta*m_eta*m_eta*Td*Td*ev-delta*delta*V(r)/np.hypot(D,2*delta*x)
        # Potential derivatives: dU/dx = x * force.
        return K+U+er+ar,K-U+ep+ap,force,(K,U,er,ar)
    def rhs(t,y):
        x,dx,z=y;r=np.exp(-t);hh=H(r-1)*HBAR
        dh=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        rho,press,force,parts=local(t,x,dx)
        return [dx,-(3+dh)*dx-x*force/(v*v*hh*hh),-3*np.exp(3*t)*press/X]
    t0=-np.log(1e16);ini=[seed,0.,0.]
    q0=np.exp(3*t0)*local(t0,seed,0.)[0]/X
    t=np.r_[np.linspace(t0,-np.log(101),300),np.linspace(-np.log(101),0,4000)[1:]]
    sol=solve_ivp(rhs,(t0,0.),ini,t_eval=t,method='DOP853',rtol=rtol,atol=[rtol*1e-10,rtol*1e-10,rtol*1e-4],max_step=.04)
    assert sol.success,sol.message
    E=[];parts=[]
    for tt,x,dx in zip(t,sol.y[0],sol.y[1]):
        rr,pp,ff,pparts=local(tt,x,dx)
        E.append(np.exp(3*tt)*rr/X);parts.append(pparts)
    E=np.asarray(E);err=E-q0-sol.y[2]
    # Compute scale using actual comoving energy excursion, not background
    # rest-mass or radiation which can hide poor conservation.
    return dict(scheme=scheme,rtol=rtol,today_x=float(sol.y[0,-1]),nfev=sol.nfev,
       max_comoving_energy_work_residual_over_X=float(max(abs(err))),
       max_relative_to_max_comoving_energy_excursion=float(max(abs(err))/max(abs(E))),
       final_energy_parts_eV4=dict(zip(['kinetic_chi','potential_chi','eta_mass_energy_increment','A_mass_energy_increment'],map(float,parts[-1]))),
       today_total_increment_fraction=float(sum(parts[-1])/rho_tot(0.)),
       scope='Fixed H, first order tiny A mass shift, independently integrated comoving work. No perturbation or thermal-origin validation.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--output',required=True);x=a.parse_args()
    row=json.loads(Path(x.input).read_text())[0];p0=scales(.1,z_on=5.,dneff=.1,muH=100.)
    p=scales(.1,z_on=5.,dneff=.1,muH=100.,X=1.05*p0['X_loop_floor'])
    out=[run(row,p,s,r) for s in ['origin-vev','vev-curvature'] for r in [2e-7,2e-9]]
    Path(x.output).write_text(json.dumps(dict(params=p,results=out),indent=2,allow_nan=False)+'\n')
    print(out,flush=True)
