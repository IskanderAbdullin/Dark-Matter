"""Frozen-PSD NR collision bracket for a narrow velocity resonance.

sigma/m = A delta(v-vr), A fixed by ONE dilute halo collision per particle
in 10 Gyr (Maxwell relative velocities). Valid for angle-independent
cross section, negligible resonance lifetime, and cold on-shell modes.
All four momenta below cut, hence a symmetric positive subdomain for
the entropy bracket. Omits rest of spectrum and actual time evolution.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,gamma,zeta
from scipy.integrate import trapezoid
from quartic_consistency import shell_nodes
from quantum_elastic_moments import lookup
from velocity_only_exposure import RHO0,C,SEC_GYR,halo_pdf
from late_activation_overlap import H

def entropy(f):
    f=np.asarray(f)
    g=np.maximum(f,1e-300)
    return np.log1p(g)+g*np.log1p(1/g)

def clock_angle(S,wlo,whi,vr,v0,rmax,A,order=8):
    """Integrate the cosmic clock over the exact NR angular interval dP."""
    # P^2+w^2=2S; both P and w share same allowed interval.
    lo=np.maximum(wlo,vr/(v0*rmax));hi=np.minimum(whi,vr/v0)
    vmax=np.sqrt(2*S)
    aa=np.arcsin(np.clip(lo/np.maximum(vmax,1e-300),1e-300,1.))
    bb=np.arcsin(np.clip(hi/np.maximum(vmax,1e-300),1e-300,1.))
    valid=(hi>lo)&(vmax>0)
    # log alpha map: dP=-w d alpha cancels the 1/w in the delta clock.
    low=np.where(valid,aa,1.);high=np.where(valid,bb,1.)
    length=np.log(high/low)
    u,w=roots_legendre(order)
    out=np.zeros_like(whi)
    for xi,wi in zip(u,w):
        a=low*np.exp(length*(xi+1)/2)
        ww=vmax*np.sin(a);r=np.clip(vr/(v0*np.maximum(ww,1e-300)),1.,rmax)
        out+=wi*RHO0*A*r**3/H(r-1)*a*length/2
    return np.where(valid,out,0.)

def calc(row,vr_kms,rho=6.700740771077128e-24,vT=62.3668181006535e5,cut=6.,order=2,angle_order=8,block=64):
    t=time.monotonic();h=row['final'];p0=h['peak_bin_q']
    e0=np.asarray(row['frozen_cell_spectrum']['edges']);f0=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    e=(e0/p0)**2
    cuts=np.array([2.,4.,cut])**2
    ee=np.unique(np.r_[0,e,cuts]);ee=ee[ee<=cut**2]
    mid=(ee[1:]+ee[:-1])/2
    Fcell=lookup(e,f0,mid)
    u,w=roots_legendre(order)
    y=mid[:,None]+np.diff(ee)[:,None]*u/2
    wt=np.diff(ee)[:,None]*w/2
    F=np.broadcast_to(Fcell[:,None],y.shape)
    # Include all incoming nodes, including empty bins, for symmetric entropy
    # quadrature; forward moments have f1*f2=0 on empty bins.
    sel=F>0
    K,wt,F=y[sel],wt[sel],F[sel]
    q,n=shell_nodes(e0,f0,12)
    N0=4*np.pi**2*n.sum()/p0**3
    K0=4*np.pi**2*(n@(q/p0)**2)/p0**3
    coldN=float(np.sum(wt*F*np.sqrt(K)))
    S0=4*np.pi**2*shell_nodes(e0,entropy(f0),12)[1].sum()/p0**3
    Tsat=(K0/(gamma(2.5)*zeta(2.5,1)))**.4
    Nsat=gamma(1.5)*zeta(1.5,1)*Tsat**1.5
    Smax=(5/3)*K0/Tsat
    assert Nsat<N0,'saturated ideal Bose upper entropy bound unavailable'
    coldK=float(np.sum(wt*F*K**1.5))
    coldS=float(np.sum(wt*entropy(F)*np.sqrt(K)))
    tc=(coldK/(gamma(2.5)*zeta(2.5,1)))**.4
    coldSmax=(5/3)*coldK/tc
    assert gamma(1.5)*zeta(1.5,1)*tc**1.5<coldN
    M=np.r_[np.sum(wt*F*np.sqrt(K)*K*K),np.sum(wt*F*np.sqrt(K)*(K>cuts[0])),np.sum(wt*F*np.sqrt(K)*(K>cuts[1]))]
    v0=p0*2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']/86*C
    rmax=min(1e7,.03*C/(cut*v0));vr=vr_kms*1e5
    A=1/(rho*10*SEC_GYR*vr*halo_pdf(vr,vT))
    i,j=np.triu_indices(len(K))
    net=np.zeros(3);absnet=np.zeros(3);cl=np.zeros(3);ent=0.;ent_abs=0.;dilute=0.
    for kk in range(0,len(i),block):
        a,b=i[kk:kk+block],j[kk:kk+block]
        S=K[a]+K[b];lo=np.maximum(0.,S-cut**2);hi=S/2
        pref=(wt[a]*wt[b]*F[a]*F[b]*np.where(a==b,1.,2.))[:,None]
        E=np.concatenate([lo[:,None],hi[:,None],np.broadcast_to(ee,(len(S),len(ee))),S[:,None]-ee[None,:]],axis=1)
        E=np.sort(np.clip(E,lo[:,None],hi[:,None]),axis=1)
        L,B=E[:,:-1],E[:,1:];d=B-L;mid=(L+B)/2
        y1,y2=K[a,None],K[b,None];p1,p2=np.sqrt(y1),np.sqrt(y2)
        f1,f2=F[a,None],F[b,None]
        for xi,wi in zip(u,w):
            y3=mid+d*xi/2;y4=S[:,None]-y3
            p3,p4=np.sqrt(np.maximum(y3,0.)),np.sqrt(np.maximum(y4,0.))
            wlo=np.maximum(abs(p1-p2),abs(p3-p4));whi=np.minimum(p1+p2,p3+p4)
            W=clock_angle(S[:,None],wlo,whi,vr,v0,rmax,A,angle_order)
            base=pref*d*wi*W
            f3,f4=lookup(e,f0,y3),lookup(e,f0,y4)
            Bfac=1+f3+f4
            ds=[y3*y3+y4*y4-y1*y1-y2*y2]
            for c in cuts[:2]:ds.append((y3>c).astype(float)+(y4>c).astype(float)-(y1>c).astype(float)-(y2>c).astype(float))
            for ix,dy in enumerate(ds):
                net[ix]+=np.sum(base*Bfac*dy)
                absnet[ix]+=np.sum(base*Bfac*abs(dy))
                cl[ix]+=np.sum(base*dy)
            dilute+=float(base.sum())
            # Positive symmetric entropy integrand on all f_i>=100.
            ok=(f1>=100)&(f2>=100)&(f3>=100)&(f4>=100)&(base>0)
            a1=np.broadcast_to(f1,ok.shape)[ok];a2=np.broadcast_to(f2,ok.shape)[ok]
            a3=f3[ok];a4=f4[ok]
            rin=np.log1p(1/a1)+np.log1p(1/a2);rout=np.log1p(1/a3)+np.log1p(1/a4)
            delta=rout-rin
            # base already includes f1*f2. Remaining product f3*f4.
            entropy_flux=.5*base[ok]*a3*a4*np.exp(np.maximum(rin,rout))*(-np.expm1(-abs(delta)))*abs(delta)
            ent+=float(entropy_flux.sum())
    return dict(vr_kms=vr_kms,halo_density=rho,halo_vT_kms=vT/1e5,halo_age_Gyr=10.,area_A_cm2_g_times_cm_s=A,
        cut_over_peak=cut,order=order,angle_order=angle_order,zmax=rmax-1,number_full=N0,cold_number_fraction=coldN/N0,
        cold_entropy_initial_per_total_particle=coldS/N0,cold_entropy_max_per_total_particle=coldSmax/N0,
        entropy_lower_over_cold_available_entropy=(ent/(4*N0))/(coldSmax-coldS),
        entropy_initial_per_particle=S0/N0,entropy_max_per_particle=Smax/N0,equilibrium_excited_fraction=Nsat/N0,
        relative_frozen_moment_changes=(net/(4*N0*M)).tolist(),dilute_relative_frozen_changes=(cl/(4*N0*M)).tolist(),
        moment_net_over_absolute=(net/np.maximum(absnet,1e-300)).tolist(),
        entropy_lower_contribution_per_particle=ent/(4*N0*N0),
        entropy_lower_over_available_entropy=(ent/(4*N0))/(Smax-S0),
        dilute_collision_exposure_cold_all_four=dilute/(2*N0*N0),
        labels=['cold_p4','N(2peak<q<cut)','N(4peak<q<cut)'],wall_seconds=time.monotonic()-t,
        scope='Frozen collisionless PSD; cold symmetric subdomain, positive entropy term only all f>=100; no actual evolution, medium resonance dressing or exclusion.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--velocities',nargs='+',type=float,default=[40.,90.,140.]);p.add_argument('--order',type=int,default=2)
    p.add_argument('--rho',type=float,default=6.700740771077128e-24);p.add_argument('--vT',type=float,default=62.3668181006535e5);
    p.add_argument('--angle-order',type=int,default=8);p.add_argument('--cut',type=float,default=6.);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];results=[]
    for v in a.velocities:
        r=calc(row,v,rho=a.rho,vT=a.vT,cut=a.cut,order=a.order,angle_order=a.angle_order)
        results.append(r);Path(a.output).write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
        print(json.dumps(r),flush=True)
