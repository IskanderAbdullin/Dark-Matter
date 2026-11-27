"""Frozen isotropic charge screening and fluctuation moments."""
import gzip
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit,zeta

C=2/np.pi**2  # four particle/antiparticle spin states per Dirac species


def raw(edges,f,m,xcut):
    p=np.asarray(edges)
    e=np.hypot(p,m)
    f=np.asarray(f)
    assert min(f)>=0. and max(f)<=1.
    # Static screening derivative split, exact for the piecewise declared cells.
    stat=C*np.diff(p*e)@f
    n=C*np.diff(p**3)@f/3
    energy=C*sum(o*quad(lambda q:q*q*np.hypot(q,m),l,h,
                        epsabs=1e-13)[0] for l,h,o in zip(p[:-1],p[1:],f))
    lim=m*xcut/np.sqrt(1-xcut*xcut)
    fp=np.r_[0.,f,0.]
    ret=C*np.sum(e[p>lim]**2*(-np.diff(fp)[p>lim]))
    h=np.maximum(p,lim)
    en=np.hypot(h,m)
    fl=C*np.diff(en**3)@(f*(1-f))/3
    return dict(static_mD2_over_g2=stat,n=n,rho=energy,noise_weight=fl,
                retarded_cut_weight=ret)


def fd(m,N,T,xcut):
    mu=brentq(lambda z:C*quad(lambda p:p*p*expit(z-np.hypot(m,p)/T),
                              0.,100.,epsabs=1e-12)[0]-N,-40.,25.)
    def noise(p):
        a=mu-np.hypot(p,m)/T
        return expit(a)*expit(-a)
    lim=m*xcut/np.sqrt(1-xcut*xcut)
    stat=C/T*quad(lambda p:p*p*noise(p),0.,100.,epsabs=1e-12)[0]
    fl=C*quad(lambda p:p*np.hypot(p,m)*noise(p),lim,100.,epsabs=1e-12)[0]
    R=C*quad(lambda p:p*p*np.hypot(p,m)*expit(mu-np.hypot(p,m)/T),
             0.,100.,epsabs=1e-12)[0]
    return dict(static_mD2_over_g2=stat,n=N,T=T,alpha=mu,rho=R,
                noise_weight=fl,retarded_cut_weight=fl/T)


def combining(P,B):
    ret=P['retarded_cut_weight']+B['retarded_cut_weight']
    S=P['static_mD2_over_g2']+B['static_mD2_over_g2']
    return dict(parent=P,daughter=B,static_mD2_vs_two_thermal_dirac=S/(2/3),
                daughter_static_fraction=B['static_mD2_over_g2']/S,
                fluctuation_over_response_effective_temperature=(P['noise_weight']+B['noise_weight'])/ret)


def run():
    check=fd(0.,3*zeta(3)/np.pi**2,1.,.7)
    # The approximate FD normalization supplied above is from 3 zeta(3)/pi^2.
    assert abs(check['static_mD2_over_g2']/(1/3)-1)<1e-7
    d=json.load(gzip.open('results/early-pole-spectrum-third-m06-2026-11-26.json.gz','rt'))
    fine=json.load(gzip.open('results/early-pole-spectrum-fine-m06-2026-11-26.json.gz','rt'))
    assert fine['x_end']==d['x_end'] and fine['mass_ratio']==d['mass_ratio']
    snaps=json.load(gzip.open('results/early-pole-separate-fermion-heat-2026-11-26.json.gz','rt'))
    x=d['x_end']; rb=d['mass_ratio']
    rows=[]
    for vel in (0.,.5,.9,.99,.999):
        P,B=(raw(d[f'endpoint_{k}_momentum_edges'],d[f'endpoint_{k}_occupations'],m,vel)
             for k,m in [('parent',x),('fermion',rb*x)])
        tb,tp=snaps['separate-fermion-v2-fast']['reference_initial_temperatures']
        rows.append(dict(case='initial_resolved',phase_velocity=vel,**combining(P,B)))
        Pf,Bf=(raw(fine[f'endpoint_{k}_momentum_edges'],fine[f'endpoint_{k}_occupations'],m,vel)
               for k,m in [('parent',x),('fermion',rb*x)])
        rows.append(dict(case='initial_resolved_other_grid',phase_velocity=vel,**combining(Pf,Bf)))
        rows.append(dict(case='initial_within_species_FD',phase_velocity=vel,
                         **combining(fd(x,P['n'],tp,vel),fd(rb*x,B['n'],tb,vel))))
        for case,s in snaps.items():
            xe=s['xend']
            rows.append(dict(case=case,phase_velocity=vel,
                         **combining(fd(xe,s['final_parent_count_over_T3'],s['final_parent_T'],vel),
                                    fd(rb*xe,s['final_daughter_count_over_T3'],s['final_daughter_T'],vel))))
    return dict(scope='Frozen isotropic collisionless charge response and independent-particle noise; '
                      'no self-consistent HTL trajectory, photon contribution to Abelian charge zero.',
                thermal_control=check,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
