"""Tree-level 3L -> (L or H)+S for the subthreshold mediator EFT.

Potential (g/2) S B^2 + lambda0 B^4/24, B=-sqrt(u)L+sqrt(1-u)H.
All four external light/heavy scalar legs attach to B vertices. Include
4 cubic-quartic diagrams and 12 three-cubic diagrams, with both internal
mass eigenstates summed coherently. The stable S lies below 2m_L.
Compute initial source only (S and H empty), no inverse evolution,
in-medium propagator, higher diagrams or bound-state matching.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid,trapezoid
from scipy.stats import qmc
from scipy.special import zeta
from soft_mixing_activation import frac_B
from direct_activation_consistency import CGS_RHO_TO_EV4,T0
from resonant_mediator_inverse import M,HBAR,C,SEC_GYR
from late_activation_overlap import H,critical_density
from velocity_only_exposure import RHO0
from two_state_halo_equilibrium import pl

def norm2(p):return p[...,0]**2-np.sum(p[...,1:]**2,axis=-1)
def stable_masses(D,delta,x,a,vT):
    # A vacuum mass M; at full mixing x=1 define S below threshold.
    shift=lambda x:-2*delta*delta*x*x/(np.hypot(D,2*delta*x)+D)
    L2=M*M+shift(x);H2=L2+np.hypot(D,2*delta*x)
    S2=4*(M*M+shift(1.))-a*M*M*(vT/C)**2
    return np.sqrt(L2),np.sqrt(H2),np.sqrt(S2),frac_B(x,delta,D)

def amplitude(p,S,u,g,lam,mL,mH,mS,outH=False):
    # p list of 4 B-like external four-vectors, S outgoing convention
    # already converted to all-incoming. Symmetric tree sum.
    out=np.zeros_like(p[0][...,0])
    for i in range(4):
        k2=norm2(S+p[i])
        GB=u/(k2-mL*mL)+(1-u)/(k2-mH*mH)
        ix=[j for j in range(4) if j!=i]
        gs=0.
        for a in range(3):
            for b in range(a+1,3):
                gs=gs+1/(norm2(p[ix[a]]+p[ix[b]])-mS*mS)
        out+=GB*(lam+g*g*gs)
    ext=u**1.5*np.sqrt(1-u if outH else u)
    return -g*ext*out

def final_state(p1,p2,p3,mL,mF,mS,cost,phi):
    # Exact CM boost, stable Kallen factor, arbitrary incoming arrays.
    E1,E2,E3=[np.sqrt(mL*mL+np.sum(p*p,axis=-1)) for p in [p1,p2,p3]]
    Et=E1+E2+E3;P=p1+p2+p3
    ss=Et*Et-np.sum(P*P,axis=-1);root=np.sqrt(ss)
    allowed=root>=mF+mS
    k2=np.maximum((ss-(mF+mS)**2)*(ss-(mF-mS)**2)/(4*ss),0.)
    k=np.sqrt(k2)
    sin=np.sqrt(np.maximum(1-cost*cost,0.))
    n=np.stack([sin*np.cos(phi),sin*np.sin(phi),cost],axis=-1)
    pcm=n*k[...,None]
    gamma=Et/root;beta=P/Et[...,None];bdot=np.sum(beta*pcm,axis=-1)
    Ef=(ss+mF*mF-mS*mS)/(2*root)
    # (gamma-1)/beta^2=gamma^2/(gamma+1).
    out=pcm+((gamma*gamma/(gamma+1)*bdot+gamma*Ef)[...,None])*beta
    eout=gamma*(Ef+bdot);PS=P-out;ES=Et-eout
    four=lambda E,p:np.concatenate([E[...,None],p],axis=-1)
    return [four(E1,p1),four(E2,p2),four(E3,p3),-four(eout,out)],-four(ES,PS),k,root,allowed,(E1,E2,E3),out,PS

def zero_rate(u,g,lam,mL,mH,mS,outH=False):
    z=np.zeros((1,3));mF=mH if outH else mL
    p,S,k,root,ok,Es,po,ps=final_state(z,z,z,mL,mF,mS,np.ones(1),np.zeros(1))
    if not ok[0]:return dict(allowed=False,coefficient=0.)
    amp=float(amplitude(p,S,u,g,lam,mL,mH,mS,outH)[0])
    coeff=amp*amp*k[0]/(4*np.pi*root[0])/(6*(2*mL)**3)
    return dict(allowed=True,coefficient=float(coeff),amplitude_eV_inv=amp,
       recoil_momentum_eV=float(k[0]),energy_release_eV=3*mL-mF-mS)

def radial_bose(alpha):
    x=np.unique(np.r_[0.,np.geomspace(1e-10,.1,16000),np.linspace(.1,10,36000)])
    f=x*x/np.expm1(x*x+alpha);cdf=cumulative_trapezoid(f,x,initial=0.)
    norm=cdf[-1];cdf/=norm
    norm_ref=np.sqrt(np.pi)/4*float(pl(1.5,alpha))
    mean=trapezoid(f*x*x,x)/norm
    mean_ref=1.5*float(pl(2.5,alpha))/float(pl(1.5,alpha))
    return x,cdf,dict(number_rel_error=norm/norm_ref-1,mean_q2_rel_error=mean/mean_ref-1)

def radial_vectors(samples,q1,q2,q3,scale):
    c2=2*samples[:,3]-1;c3=2*samples[:,4]-1;phi=2*np.pi*samples[:,5]
    si2=np.sqrt(1-c2*c2);si3=np.sqrt(1-c3*c3)
    z=np.zeros_like(q1)
    p1=scale*np.stack([z,z,q1],axis=-1)
    p2=scale*np.stack([q2*si2,z,q2*c2],axis=-1)
    p3=scale*np.stack([q3*si3*np.cos(phi),q3*si3*np.sin(phi),q3*c3],axis=-1)
    return p1,p2,p3,2*samples[:,6]-1,2*np.pi*samples[:,7]

def avg_rate(samples,radial,scale,u,g,lam,mL,mH,mS,final_f,outH=False):
    q1,q2,q3=[radial(samples[:,i]) for i in range(3)]
    p1,p2,p3,c,phi=radial_vectors(samples,q1,q2,q3,scale)
    mF=mH if outH else mL
    p,S,k,root,ok,Es,po,ps=final_state(p1,p2,p3,mL,mF,mS,c,phi)
    if not np.any(ok):return dict(coefficient=0.,classical_coefficient=0.,allowed_fraction=0.,max_amplitude=0.)
    # Exclude closed channels before evaluating propagators.
    p=[v[ok] for v in p];S=S[ok];Es=[v[ok] for v in Es]
    amp=amplitude(p,S,u,g,lam,mL,mH,mS,outH)
    base=amp*amp*k[ok]/(4*np.pi*root[ok])/(48*Es[0]*Es[1]*Es[2])
    boost=np.ones_like(base) if outH else 1+final_f(np.linalg.norm(po[ok],axis=-1))
    W=base*boost
    return dict(coefficient=float(W.sum()/len(samples)),classical_coefficient=float(base.sum()/len(samples)),
       allowed_fraction=float(ok.mean()),max_amplitude=float(max(abs(amp))),
       weighted_final_S_p_eV=float(W@np.linalg.norm(ps[ok],axis=-1)/max(W.sum(),1e-300)),
       weighted_Bose_factor=float(W.sum()/max(base.sum(),1e-300)))

def halo(p,V=44.1,alpha=1e-6,powers=(16,18),replicates=4):
    vT=np.sqrt(2)*V*1e5;scale=M*vT/C;rho=critical_density(vT);n=rho*CGS_RHO_TO_EV4/M
    u=p['u'];D=p['D'];delta=p['delta'];a=p['a'];g=p['g_eV']
    mL,mH,mS,u=stable_masses(D,delta,1.,a,vT);lam=3*g*g/(mS*mS)
    x,cdf,meta=radial_bose(alpha)
    radial=lambda v:np.interp(v,cdf,x)
    final=lambda q:np.exp(-alpha-(q/scale)**2)/(-np.expm1(-alpha-(q/scale)**2))
    out=[]
    for power in powers:
        for rep in range(replicates):
            s=qmc.Sobol(d=8,scramble=True,seed=137+rep).random_base2(power)
            L=avg_rate(s,radial,scale,u,g,lam,mL,mH,mS,final,False)
            Hh=avg_rate(s,radial,scale,u,g,lam,mL,mH,mS,final,True)
            coeff=L['coefficient']+Hh['coefficient']
            # L+S final removes 2 L; H+S removes 3 L. Count real S events separately.
            out.append(dict(power=power,rep=rep,L=L,H=Hh,initial_S_events_per_particle_10Gyr=coeff*n*n*10*SEC_GYR/HBAR,
              initial_L_depletion_10Gyr=(2*L['coefficient']+3*Hh['coefficient'])*n*n*10*SEC_GYR/HBAR))
    return dict(params=p,V0_kms=V,masses_eV=[mL,mH,mS],u=u,lambda0=lam,alpha=alpha,
       radial_checks=meta,zero_L=zero_rate(u,g,lam,mL,mH,mS),zero_H=zero_rate(u,g,lam,mL,mH,mS,True),
       runs=out,scope='Initial dilute S,H source with initial BE L including final L stimulation. 3->2 tree diagrams; no inverse, finite density mediator dressing or nonlinear evolution.')

def cosmic_zero(row,p,rr):
    r=rr['r'];x=rr['x'];dx=rr['dx'];t=rr['t'];keep=r<=101
    r,x,dx,t=r[keep],x[keep],dx[keep],t[keep]
    vT=np.sqrt(2)*44.1e5;D=p['D'];delta=p['delta'];a=p['a'];g=p['g_eV']
    tr=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']
    e=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    Cc=[];mom=[];beta=[];rate=[];rH=[]
    for R,xx in zip(r,x):
        ml,mh,ms,u=stable_masses(D,delta,xx,a,vT);lam=3*g*g/ms**2
        L=zero_rate(u,g,lam,ml,mh,ms);Hh=zero_rate(u,g,lam,ml,mh,ms,True)
        n=RHO0*CGS_RHO_TO_EV4/M*R**3
        hh=HBAR*H(R-1)
        q=L['recoil_momentum_eV']/(tr*R)
        j=np.searchsorted(e,q)-1;ff=f[j] if 0<=j<len(f) else 0.
        coeff=L['coefficient']*(1+ff)
        Cc.append(coeff);mom.append(L['recoil_momentum_eV']);beta.append(ff)
        rate.append(coeff*n**3/hh);rH.append(Hh['coefficient']*n*n/hh)
    mom=np.array(mom);rate=np.array(rate);n=RHO0*CGS_RHO_TO_EV4/M*r**3
    sweep=1+np.gradient(np.log(mom),t)
    fS0=2*np.pi**2*rate/mom**3/np.maximum(abs(sweep),1e-300)
    return dict(initial_S_fraction=float(trapezoid(rate/n,t)),initial_HS_event_fraction=float(trapezoid(rH,t)),
        maximum_zero_width_line_occupation=float(max(fS0)),minimum_abs_line_sweep=float(min(abs(sweep))),
        maximum_final_L_occupation=float(max(beta)),max_recoil_over_peak_today=float(max(mom/(tr*r*row['final']['peak_bin_q']))),
        scope='Cold zero-incoming-momentum initial-source approximation. S occupation from swept monochromatic line; no inverse/stimulation, hot-tail or stationary-line regularization.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--tradeoff',required=True);ap.add_argument('--output',required=True);ap.add_argument('--powers',nargs='+',type=int,default=[16,18]);ap.add_argument('--replicates',type=int,default=3);a=ap.parse_args()
    data=json.loads(Path(a.tradeoff).read_text());out=[]
    for b in data['exact']:
        if b['z_on'] not in [8.,20.]:continue
        p=dict(u=b['u'],D=b['D_eV2'],delta=b['delta_eV2'],a=b['a'],g_eV=b['g_eV'],z_on=b['z_on'])
        r=halo(p,powers=a.powers,replicates=a.replicates);out.append(r)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print(p,r['zero_L'],r['zero_H'],[(x['power'],x['rep'],x['initial_L_depletion_10Gyr'],x['L']['weighted_Bose_factor'],x['H']['allowed_fraction']) for x in r['runs']],flush=True)
