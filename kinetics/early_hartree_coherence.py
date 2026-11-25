"""Coherent post-portal flavor evolution with self-consistent repulsive B Hartree.

Prescribed homogeneous early trigger x=seed, two nearly degenerate scalars,
common positive energy E=sqrt(p^2+m^2), H_flavor=(mass^2-common)/(2E).
Evolve normalized complex amplitudes for each occupied comoving shell.
Small mixing is scaled out of the heavy amplitude for accurate integration.
No diagonal scattering or dephasing in this pilot; evaluate collision and
common-energy validity diagnostics along the trajectory. This is not a full
finite-temperature quantum kinetic equation.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,trapezoid
from scipy.optimize import brentq
from quartic_consistency import shell_nodes
from soft_mixing_activation import scales,frac_B
from late_activation_overlap import H
from direct_activation_consistency import CGS_RHO_TO_EV4,T0
from resonant_mediator_inverse import M,HBAR
from velocity_only_exposure import RHO0

N0=RHO0*CGS_RHO_TO_EV4/M

def samples(row,order=2,trim=0.):
    e=np.asarray(row['frozen_cell_spectrum']['edges'])
    f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    q,w=shell_nodes(e,f,order);f=np.repeat(f,order)
    w=w/w.sum()
    # Optional remove only an explicitly reported small low-momentum tail,
    # chosen using both number and ultrarelativistic variance (w/q).
    ix=np.argsort(q);q,w,f=q[ix],w[ix],f[ix]
    active=w>0;q,w,f=q[active],w[active],f[active]
    ncut=np.cumsum(w);vcut=np.cumsum(w/q)/sum(w/q)
    keep=(ncut>=trim)|(vcut>=trim)
    info=dict(dropped_number=float(sum(w[~keep])),dropped_UR_variance=float(sum(w[~keep]/q[~keep])/sum(w/q)))
    return q[keep],w[keep],f[keep],info

def equilibrium(E,w,r,D,off,lam):
    # Pure lower eigenstate, B variance = n0*r^3*<u_eff/E>.
    V=N0*r**3*np.sum(w/E)
    def f(Deff):
        u=frac_B(1.,off,Deff)
        return Deff-D-lam*V*u/2
    Deff=brentq(f,D,D+lam*V/2+1e-20,xtol=1e-17,rtol=1e-13)
    return Deff,frac_B(1.,off,Deff)

def run(row,p,seed=1e-8,rstart=1e14,rend=3e11,order=2,trim=1e-10,rtol=1e-8,initial='A',points=12000):
    start=time.monotonic()
    q,w,f,meta=samples(row,order,trim)
    Tr=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']
    D=p['Delta'];off=p['delta']*seed;lam=p['lambda_B'];eps=off/D
    n=len(q);t0=-np.log(rstart);t1=-np.log(rend)
    r=rstart;E=np.hypot(M,Tr*q*r);Deq,ueq=equilibrium(E,w,r,D,off,lam)
    if initial=='A':b=np.zeros(n,complex)
    elif initial=='lower':b=np.full(n,-np.sqrt(ueq),complex)
    elif initial=='phase+':b=np.full(n,1j*eps,complex)
    elif initial=='phase-':b=np.full(n,-1j*eps,complex)
    else:raise ValueError(initial)
    a=np.sqrt(1-abs(b)**2).astype(complex);c=b/eps
    # Constant positive energy offset chosen to make the local vacuum
    # quadratic form >= D/4 * f_B. Conformal energy at fixed chi cannot
    # increase under expansion, even with nonlinear Hartree oscillations.
    shift=off*off/D
    W0=N0*rstart*np.sum(w/E)
    escale=shift*W0
    def local(t,y):
        r=np.exp(-t);E=np.hypot(M,Tr*q*r);h=HBAR*H(r-1)
        a=y[:n];b=eps*y[n:2*n]
        B=abs(b)**2;A=abs(a)**2;co=(a*np.conj(b)).real
        weight=N0*w*r/E
        JB=weight@B
        V=r*r*JB
        Deff=D+lam*V/2
        Omega=off/(2*E*h);WW=Deff/(2*E*h)
        vv=D/2*B+off*co+shift*(A+B)
        energy=weight@vv+lam*r*r*JB*JB/8
        dw=-weight*(M/E)**2
        dJB=dw@B
        work=dw@vv+lam*r*r*(-JB*JB+JB*dJB)/4
        return a,b,Omega,WW,energy,work,V,Deff,E,JB
    def rhs(t,y):
        a,b,O,W,energy,work,V,Deff,E,JB=local(t,y)
        return np.r_[-1j*O*b,-1j*(W*(b/eps)+O/eps*a),work/escale]
    y0=np.r_[a,c,0.+0.j]
    E0=local(t0,y0)[4]
    ts=np.linspace(t0,t1,points+1)
    sol=solve_ivp(rhs,(t0,t1),y0,t_eval=ts,rtol=rtol,atol=rtol*.02,method='DOP853',max_step=.03)
    assert sol.success,sol.message
    stats=[];last=0.
    maxnorm=0.;maxerr=0.;maxf=0.;maxH=0.;maxvh=0.;maxcommon=0.;maxbound=0.;maxgain=0.
    gamma=[];modegain=[];Hfrac=[];dets=[];energ=[]
    for i,t in enumerate(ts):
        y=sol.y[:,i];r=np.exp(-t)
        a,b,O,W,energy,work,V,Deff,E,JB=local(t,y)
        B=abs(b)**2;A=abs(a)**2
        u=frac_B(1.,off,Deff)
        heavy=abs(np.sqrt(u)*a+np.sqrt(1-u)*b)**2
        FH=float(max(f*B));frac=float(w@B)
        actualH=float(w@heavy)
        h=HBAR*H(r-1)
        maxnorm=max(maxnorm,float(max(abs(A+B-1))))
        err=float(abs(energy-E0-escale*y[-1].real)/escale)
        maxerr=max(maxerr,err);maxf=max(maxf,FH);maxH=max(maxH,actualH)
        maxvh=max(maxvh,float((Deff-D)/D));maxcommon=max(maxcommon,float(max((Deff+abs(off))/(E*E))))
        # Hartree energy provides a trajectory-independent coherent bound
        # on J_B while chi is constant and there are no collisions.
        Jbound=4*E0/D
        maxbound=max(maxbound,float(JB/Jbound))
        # Overestimate per-particle dilute B collision rate, including
        # two final-state occupation factors. Actual I_B is a quadrature
        # on the original shells; its use is diagnostic, not a bootstrap.
        tau=lam*lam*V*(1+FH)**2/(32*np.pi*M*h)
        gamma.append(tau)
        IB=2*np.pi**2*N0*r**3*np.sum(w*B/(E*(q*Tr*r)))
        gain=lam*lam*(1+2*FH)*IB*IB/(64*np.pi**3*M*h)
        modegain.append(gain)
        Hfrac.append(actualH);dets.append((Deff-D)/D);energ.append(energy/escale)
        if i%(max(1,points//500))==0 or i==points:
            stats.append(dict(r=float(r),B_fraction=frac,H_eigen_fraction=actualH,max_fB=FH,
              hartree_over_split=float((Deff-D)/D),u_inmedium=float(u),W_over_H_min=float(min(W)),W_over_H_peak=float(W[np.argmax(f)]),
              energy_over_initial=float(energy/E0),energy_work_residual_over_scale=err,
              B_collision_over_H_upper=float(tau),B_mode_gain_over_H_diag=float(gain)))
    r=rend;E=np.hypot(M,Tr*q*r);Deq,ueq=equilibrium(E,w,r,D,off,lam)
    a,b,O,W,energy,work,V,Deff,E,JB=local(t1,sol.y[:,-1])
    u=frac_B(1.,off,Deff)
    heavy=abs(np.sqrt(u)*a+np.sqrt(1-u)*b)**2
    # Use actual trajectory angle derivative as a local adiabaticity
    # diagnostic. Fast coherent Hartree oscillations can violate the
    # smooth-angle approximation even if their amplitude is small.
    angle=.5*np.arctan2(2*off,D*(1+np.array(dets)))
    g=np.gradient(angle,ts)
    Wpeak=(D*(1+np.asarray(dets)))/(2*np.hypot(M,Tr*q[np.argmax(f)]*np.exp(-ts))*HBAR*H(np.exp(-ts)-1))
    return dict(params=p,seed=seed,initial=initial,rstart=rstart,rend=rend,order=order,trim=trim,rtol=rtol,points=points,
       n_modes=n,meta=meta,nfev=sol.nfev,wall_seconds=time.monotonic()-start,
       max_norm_residual=maxnorm,max_conformal_work_residual_over_scale=maxerr,
       max_B_occupation=maxf,max_H_eigen_fraction=maxH,max_Hartree_over_vacuum_split=maxvh,
       max_common_energy_expansion_parameter=maxcommon,max_JB_over_energy_bound=maxbound,
       integrated_B_collision_per_particle_upper_on_trajectory=float(trapezoid(gamma,ts)),
       integrated_B_mode_gain_diag_on_trajectory=float(trapezoid(modegain,ts)),
       final_eigen_fraction=float(w@heavy),final_B_fraction=float(w@abs(b)**2),
       final_Hartree_over_split=float((Deff-D)/D),final_ground_state_Hartree_over_split=float((Deq-D)/D),
       final_angle_rate_over_peak_gap=float(abs(g[-1]/Wpeak[-1])),
       late_max_angle_rate_over_peak_gap=float(max(abs(g[ts>=ts[-1]-.25]/Wpeak[ts>=ts[-1]-.25]))),
       final=dict(q=q.tolist(),weights=w.tolist(),f_initial=f.tolist(),a_real=a.real.tolist(),a_imag=a.imag.tolist(),
         b_real=b.real.tolist(),b_imag=b.imag.tolist(),heavy_probability=heavy.tolist()),
       history=stats,
       scope='Common-energy two-state Hartree Hamiltonian, constant early chi, collisionless amplitudes with B self-energy; no collision-induced redistribution, thermal portal after rstart, offdiagonal dissipative terms or proof of full late adiabaticity.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--sigma',type=float,default=.3);ap.add_argument('--seed',type=float,default=1e-8);ap.add_argument('--rstart',type=float,default=1e14);ap.add_argument('--rend',type=float,default=3e11)
    ap.add_argument('--order',type=int,default=2);ap.add_argument('--trim',type=float,default=1e-10);ap.add_argument('--rtol',type=float,default=1e-8)
    ap.add_argument('--initial',choices=['A','lower','phase+','phase-'],default='A');ap.add_argument('--points',type=int,default=12000);a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0]
    p=scales(a.sigma,z_on=8 if a.sigma!=.1 else 5.,muH=30 if a.sigma==.3 else 100.,dneff=.1)
    r=run(row,p,a.seed,a.rstart,a.rend,a.order,a.trim,a.rtol,a.initial,a.points)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k not in ['history','final','params']}),flush=True)
