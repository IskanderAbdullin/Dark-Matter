"""Single-channel narrow elastic Breit-Wigner realization and inverse-decay test.

Identical scalar s wave: peak sigma=8*pi/k^2, E=m v_rel^2/(4c^2).
For a pure elastic resonance the required halo area fixes its tiny vacuum
width, hence tests whether an instantaneous elastic kernel is appropriate.
The explicit R occupation equation uses the exact nonrelativistic decay
angle integral for fixed A PSD, includes final-state Bose stimulation and
inverse reactions, and stops interpreting beyond significant backreaction.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from quartic_consistency import shell_nodes
from velocity_only_exposure import halo_pdf,C,RHO0,SEC_GYR
from late_activation_overlap import H

M=86.
HBAR=6.582119569e-16
HBARC=1.973269804e-5
MASSG=M*1.7826619216278976e-33

def width(vr_kms,rho=6.700740771077128e-24,vT=62.3668181006535e5,age=10.,events=1.,branch=1.):
    vr=vr_kms*1e5
    A=events/(rho*age*SEC_GYR*vr*halo_pdf(vr,vT))
    k=M*vr/(2*C)
    E=k*k/M
    peak=8*np.pi*HBARC**2/k**2/MASSG
    dedv=M*vr/(2*C*C)
    G=2*A*dedv/(np.pi*peak*branch**2)
    gel=G*branch
    beta=vr/(2*C)
    cubic=np.sqrt(gel*32*np.pi*(2*M+E)/beta)
    life=HBAR/G
    x=age*SEC_GYR/life
    completion=1+np.expm1(-x)/x
    return dict(vr_kms=vr_kms,rho=rho,vT_kms=vT/1e5,age_Gyr=age,halo_events=events,
        area=A,peak_unitarity_sigma_over_m=peak,excess_mass_eV=E,width_eV=G,elastic_width_eV=gel,
        width_over_excess=G/E,velocity_HWHM_cm_s=G/(2*dedv),lifetime_Gyr=life/SEC_GYR,
        uniform_halo_production_completion_fraction=completion,cubic_g_eV=cubic,elastic_branching=branch)

def averages(Q,b,e,f):
    # Exact piecewise-constant angular average over c=|cos theta| in [0,1].
    a=Q*Q/4+b*b;d=Q*b
    v=np.maximum(d[:,None],1e-300)
    edges=np.sort(np.clip(np.c_[np.zeros(len(Q)),np.ones(len(Q)),(e[None,:]-a[:,None])/v,
                           (a[:,None]-e[None,:])/v],0.,1.),axis=1)
    dw=np.diff(edges,axis=1);mid=(edges[:,1:]+edges[:,:-1])/2
    def occ(y):
        j=np.searchsorted(e,y,side='right')-1
        return np.where((j>=0)&(j<len(f)),f[np.clip(j,0,len(f)-1)],0.)
    f1=occ(a[:,None]+v*mid);f2=occ(a[:,None]-v*mid)
    return np.sum(dw*f1*f2,axis=1),1+np.sum(dw*(f1+f2),axis=1)

def frozen_inverse(row,params,cut=6.,nQ=150,nz=600):
    start=time.monotonic();p0=row['final']['peak_bin_q'];theta=row['final']['theta']
    e0=np.asarray(row['frozen_cell_spectrum']['edges']);f0=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    N=float(shell_nodes(e0,f0,12)[1].sum())/p0**3
    e=(e0/p0)**2
    ee=np.unique(np.r_[0.,e[e<cut**2],cut**2]);mid=(ee[:-1]+ee[1:])/2
    j=np.searchsorted(e,mid,side='right')-1
    f=np.where((j>=0)&(j<len(f0)),f0[np.clip(j,0,len(f0)-1)],0.)
    v0=p0*2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/theta/M*C
    vr=params['vr_kms']*1e5
    rmax=min(1e7,.03*C/(cut*v0));rmin=max(1.,vr/(2*v0*cut))
    edges=np.r_[0.,np.geomspace(1e-4,2*cut,nQ+1)]
    Q,w=shell_nodes(edges,np.ones(len(edges)-1),2)
    y=np.zeros_like(Q);l=np.linspace(-np.log(rmax),-np.log(rmin),nz+1)
    gamma=params['elastic_width_eV']/HBAR
    hist=[];first={};max_stim=0.
    for i,(l0,l1) in enumerate(zip(l[:-1],l[1:])):
        r=np.exp(-(l0+l1)/2);b=vr/(2*v0*r)
        F,G=averages(Q,b,ee,f)
        dt=(l1-l0)/H(r-1)
        change=-np.expm1(-gamma*G*dt)
        y+=(F/G-y)*change
        frac=2*float(w@y)/N
        stim=float(np.max(gamma*G/H(r-1)))
        max_stim=max(max_stim,stim)
        for lim in (.001,.01,.1,1.):
            if frac>=lim and str(lim) not in first:
                first[str(lim)]=dict(z=float(np.exp(-l1)-1),fraction=frac,maximum_stimulated_width_over_H=stim)
        if i%10==0 or i==len(l)-2:
            hist.append(dict(z=float(np.exp(-l1)-1),fraction_if_A_frozen=frac,
               max_fR=float(y.max()),max_stimulated_decay_over_H=stim))
    return dict(cut_over_peak=cut,nQ=nQ,nz=nz,peak_velocity_today_cm_s=v0,
        start_z=rmax-1,end_z=rmin-1,maximum_vacuum_width_over_H=float(gamma/H(rmin-1)),
        maximum_stimulated_width_over_H=max_stim,first_crossings=first,
        maximum_frozen_fraction=max(h['fraction_if_A_frozen'] for h in hist),
        final_frozen_fraction=hist[-1]['fraction_if_A_frozen'],history=hist,
        Q=Q.tolist(),fR=y.tolist(),wall_seconds=time.monotonic()-start,
        scope='R explicitly evolved but A frozen; fraction above percent requires depletion and is NOT physical abundance. Only cold A subset included; no R interactions or medium dispersion.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--velocities',nargs='+',type=float,default=[40.,90.,140.]);p.add_argument('--nQ',type=int,default=150);p.add_argument('--nz',type=int,default=600)
    p.add_argument('--rho',type=float,default=6.700740771077128e-24);p.add_argument('--vT',type=float,default=62.3668181006535e5);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];out=[]
    for vr in a.velocities:
        r=width(vr,a.rho,a.vT);r['inverse']=frozen_inverse(row,r,nQ=a.nQ,nz=a.nz)
        out.append(r);Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print(json.dumps({k:v for k,v in r.items() if k!='inverse'}),{k:v for k,v in r['inverse'].items() if k not in ('history','Q','fR')},flush=True)
