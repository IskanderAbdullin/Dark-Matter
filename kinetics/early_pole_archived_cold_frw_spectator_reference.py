"""Comoving live cold window with an independently redshifting exterior."""
import argparse
import json
from types import MethodType

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import roots_legendre
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_joint_hard_packets_frw_reference import Expanding
from early_pole_frozen_dynamic_vlasov_reference import load


def run(grid='fine',duration=8e-5,wide=1,n=8,level=0,charge_floor=.45,field_floor=.45,expanding_only=False,clock="both"):
    C=Cold(grid,wide,level,n,True,charge_floor,field_floor,duration)
    M=C.M
    M.weights=MethodType(Expanding.weights,M)
    M.eta=0.
    L=M.norm
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    spec=[]
    z,w=roots_legendre(8)
    for i,(lab,m) in enumerate([('parent',.05),('fermion',.003)],1):
        es=np.array(d[f'endpoint_{lab}_momentum_edges'])
        f=np.array(d[f'endpoint_{lab}_occupations'])
        def integral_set(lower,upper):
            lo=np.maximum(es[:-1],lower)
            hi=np.maximum(lo,np.minimum(es[1:],upper))
            half=(hi-lo)[:,None]/2
            nu=(hi+lo)[:,None]/2+half*z
            mass=half*w*nu*nu*f[:,None]*2/np.pi**2
            mask=mass.ravel()>0.
            return nu.ravel()[mask],mass.ravel()[mask]
        high=M.edge_sets[i][-1]
        low=M.edge_sets[i][0]
        def basis(floor):
            q0,c0=integral_set(floor,low)
            q1,c1=integral_set(high,np.inf)
            return np.r_[q0,q1],np.r_[c0,c1]
        spec.append((basis(field_floor),basis(0.),m))
    photon_integrand=lambda p:p**3/np.expm1(p)/np.pi**2 if p>0. else 0.
    spectator_ph=(quad(photon_integrand,0.,M.edge_sets[0][0])[0]+
                  quad(photon_integrand,M.edge_sets[0][-1],80.)[0])
    rho_vis=np.pi**2*106.75/30

    def exterior(eta):
        D=np.exp(2*eta)
        J=[M.background[0]];jd=[0.];R=spectator_ph;rd=0.
        for (p,f),(q,u),m in spec:
            E=np.sqrt(p*p+D*m*m)
            J.append((np.pi**2/2)*np.sum(f/E))
            jd.append(-(np.pi**2/2)*D*m*m*np.sum(f/E**3))
            E=np.sqrt(q*q+D*m*m)
            R+=u@E
            rd+=D*m*m*(u@np.reciprocal(E))
        return np.array(J),np.array(jd),R,rd
    M.spectators=exterior
    assert max(abs(exterior(0.)[0]-M.background))<3e-12

    Y=np.r_[M.init.ravel(),C.ns,0.,0.,0.]
    def total(y,eta):
        M.eta=eta
        _,_,_,E,S=M.weights(y[:C.j0].reshape(-1,2),True)
        _,r,s,*_=C.scalar(y[C.j0:C.j0+C.n])
        return E+r,S+s,M.trace
    tmax=duration
    R=[]
    clocks=[(0.,False),(tmax,False),(tmax,True)]
    if expanding_only:clocks=clocks[1:]
    if clock!="both":clocks=[(eta,H) for eta,H in clocks if H==(clock=="full")]
    for etaend,full_H in clocks:
        for pair in (False,True):
            e,s,tr=total(Y,0.)
            h=2e-5
            fd=(total(Y,h)[0]-total(Y,-h)[0])/(2*h)
            assert abs(fd-tr)<2e-10
            ratio=etaend/tmax
            def rhs(t,y):
                M.eta=t*ratio
                M.norm=L*np.exp(-M.eta)
                z=C.rhs(y[:-1],pair)
                if full_H:
                    rho=C.last_total_energy
                    z=z/np.sqrt(1+rho/rho_vis)
                return np.r_[z,M.trace*ratio]
            sol=solve_ivp(rhs,(0.,tmax),Y,method='DOP853',rtol=2e-10,atol=3e-13,max_step=3e-5)
            assert sol.success
            a=sol.y[:,-1]
            e1,s1,_=total(a,etaend)
            resid=e1-e-a[-1]
            assert abs(resid-a[-2])<2e-13 and abs(resid)<2e-13
            assert abs(s1-s-a[-3])<3e-11
            assert abs(np.sum((a-Y)[:C.j0].reshape(-1,2)[:,0]))<1e-13
            R.append(dict(grid=grid,end_eta=etaend,pair=pair,
                          scalar_gain=np.sum((a-Y)[C.j0:C.j0+C.n]),energy_drift=resid,
                          integrated_work=a[-2],redshift_integral=a[-1],trace_fd_error=fd-tr,
                          use_conditional_full_H=full_H, initial_dark_to_visible=e/rho_vis,
                          final_dark_to_visible=e1/rho_vis))
    return dict(grid=grid,wide=wide,level=level,n=n,charge_floor=charge_floor,field_floor=field_floor,
                exposure=duration,scope='Tagged finite cold support, frozen outside occupations that redshift, '
                      'reference vacuum weights at initial scale. No final physical relic.',rows=R)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=['fine','third'],default='fine')
    p.add_argument('--exposure',type=float,default=8e-5)
    p.add_argument('--wide',type=int,default=1,choices=[0,1,2])
    p.add_argument('--n',type=int,default=8)
    p.add_argument('--level',type=int,default=0)
    p.add_argument('--floor',type=float,default=.45)
    p.add_argument('--field-floor',type=float,default=.45)
    p.add_argument('--expanding-only',action='store_true')
    p.add_argument('--clock',choices=['both','reference','full'],default='both')
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.exposure,a.wide,a.n,a.level,a.floor,a.field_floor,a.expanding_only,a.clock),indent=2))
