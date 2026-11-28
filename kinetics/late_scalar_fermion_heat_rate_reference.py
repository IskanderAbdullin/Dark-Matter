"""Below-pole FD/BE heat conductance near conditional thermal saturation."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit

from continuum_dark_normal_capacity_reference import fugacity,moments,rule
from soft_scalar_elastic_graph_reference import spin
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M


def state_nodes(T,mass,eta,bose,n,bath_T=None):
    x,w=rule(n)
    e=np.sqrt([0.,.001,.01,.06,.2,.6,1.5,4.,8.,16.,32.,60.])
    y=(e[:-1,None]*(1-x)+e[1:,None]*(1+x))/2
    k=(T*y*y).ravel(); E=k+mass; p=np.sqrt(k*(k+2*mass))
    x=k/(T if bath_T is None else bath_T)-eta
    f=np.exp(-x)/(-np.expm1(-x)) if bose else expit(-x)
    w=(w*(e[1:]-e[:-1])[:,None]*T*y).ravel()
    w*=p*E*(1 if bose else 4)*f/(2*np.pi**2)
    return k,E,p,w


def rates(T,N,ep,ef,n=16,mu=12,out=16,az=4,Tphi=None):
    diff_beta=0. if Tphi is None else 1/Tphi-1/T
    pp=state_nodes(T,mp,ep,True,n,bath_T=Tphi)
    cc=state_nodes(T,m,ef,False,n)
    c,w=rule(mu); cf,wf=rule(out)
    phi=(np.arange(az)+.5)*2*np.pi/az
    K,E,p,kin=cc
    result=np.zeros(3)
    for k,F,u,ui in zip(*pp):
        diff=(k*E[:,None]+mp*K[:,None])-u*p[:,None]*c
        threshold=(m+mp)**2
        s=threshold+2*diff; root=np.sqrt(s)
        sq=2*diff*(4*m*mp+2*diff)/(4*s)
        q=np.sqrt(np.maximum(sq,0.))
        L=E[:,None]+F
        P=np.sqrt(np.maximum(0.,L*L-s))
        B=P*q/root
        rest=(mp*(m+mp)+diff)/root
        center=L*rest/root
        ini=np.clip(np.divide(F-center,B,out=np.zeros_like(B),where=B>0),-1.,1.)
        high=center+B
        lo=(rest*rest+(P*mp/root)**2)/high
        change=(lo-F)[:,:,None]+B[:,:,None]*(cf+1)
        Ep=F+change
        Ec=L[:,:,None]-Ep
        sout=np.sqrt(np.maximum(0.,1-cf*cf))
        tot=0.
        for azimuth in np.cos(phi):
            angle=ini[:,:,None]*cf+np.sqrt(np.maximum(0.,1-ini*ini))[:,:,None]*sout*azimuth
            tr=-2*sq[:,:,None]*(1-angle)
            # This calculation is for the very low-T branch far from the pole.
            assert np.max(s)<.82*M*M
            tot+=spin(s[:,:,None],tr)/az
        assert np.min(tot)>=0.
        zphi=np.maximum((Ep-mp)/(T if Tphi is None else Tphi),np.finfo(float).tiny)-ep
        quantum=expit((Ec-m)/T-ef)/(-np.expm1(-zphi))
        v=q*root/(F*E[:,None])
        coeff=kin[:,None]*w/2*v/(32*np.pi*s)
        W=ui*coeff[:,:,None]*wf*tot*quantum
        heat=np.sum(W*(-np.expm1(-change*diff_beta))*change)/2
        result+=np.array([W.sum(),np.sum(W*change*change),heat])
    if Tphi is not None:
        return float(result[2])
    cv=[]
    for j,(mass,bose) in enumerate(((mp,True),(m,False))):
        base=moments(T,mass,ep if bose else ef,bose)[1]
        T1=T*(1+1e-4)
        e1=fugacity(T1,mass,N[j],bose,20)
        cv.append((moments(T1,mass,e1,bose)[1]-base)/(T1-T))
    assert min(cv)>0.
    heat=result[1]/(2*T*T)*sum(1./np.array(cv))
    return dict(scalar_scattering_rate_GeV=float(result[0]/N[0]),
                two_temperature_heat_relaxation_GeV=float(heat),
                equilibrium_energy_conductance_GeV4=float(result[1]/(2*T*T)))


def compare(cap,dil):
    n=np.array(cap['species_counts_GeV3'])
    base=cap['actual_visible_T_GeV']
    crossing=dil['crossings'][0]
    nodes=[min(dil["rows"],key=lambda q:abs(np.log(q["scale_factor_relative_to_snapshot"]/b)))
           for b in (11.,22.,45.,100.,200.)]+[crossing]
    ans=[]
    for q in nodes:
        b=q['scale_factor_relative_to_snapshot']
        T=q['reversible_common_T_GeV']
        N=n/b**3
        ep=np.log(q['reversible_scalar_fugacity'])
        ef=fugacity(T,m,N[1],False,20)
        r=rates(T,N,ep,ef,12,10,14,4)
        fine=rates(T,N,ep,ef,18,14,20,6)
        h=np.sqrt(np.pi**2/45)*(base/b)**2/MPL
        err=max(abs(r[k]/fine[k]-1.) for k in r)
        assert err<.025
        jvp=None
        if not ans:
            eps=2e-5
            t0,t1=T*(1-eps),T*(1+eps)
            x=[rates(T,N,fugacity(p,mp,N[0],True,20),ef,12,10,14,4,Tphi=p) for p in (t0,t1)]
            jvp=abs((x[0]-x[1])/(t1-t0)/r["equilibrium_energy_conductance_GeV4"]-1.)
            assert jvp<5e-7
        ans.append(dict(fixed_count_temperature_derivative_defect=jvp,b=b,T_dark_common_GeV=T,relative_quadrature_change=err,photon_only_H_GeV=h,
                        scalar_scattering_over_photon_H=fine["scalar_scattering_rate_GeV"]/h,
                        heat_relaxation_over_photon_H=fine["two_temperature_heat_relaxation_GeV"]/h))
    return ans


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("capacity")
    p.add_argument("dilution")
    a=p.parse_args()
    cap=json.loads(Path(a.capacity).read_text())['rows'][-1]
    d=json.loads(Path(a.dilution).read_text())
    print(json.dumps(dict(rows=compare(cap,d)),indent=2,allow_nan=False))
