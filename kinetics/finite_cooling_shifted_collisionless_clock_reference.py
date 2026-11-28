"""Correlated expansion bounds for the shifted, frozen conversion proxy."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from finite_cooling_collisionless_expansion_reference import Fluid
from finite_cooling_reheat_portal_reference import SourceFold, TF
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m, mp


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


class Shifted:
    def __init__(self,clock,edge,prompt,redges,pair,ref,fac):
        self.clock=clock
        self.prompt=Fluid(clock,edge,prompt,fac)
        self.fac=fac
        self.omega=float(np.sum(np.asarray(pair["new_scalar_omega"])[:,0]))
        assert abs(self.omega/ref["new_gross_early_scalar_omega"]-1.)<1e-12
        self.fraction=self.omega*m/(mp*prompt["species"][1]["Omega_h2"])
        assert 0.<self.fraction<.002
        shape=pair["new_scalar_recoil_shape"]
        self.z=redges**2
        w=np.sum(shape["count_fraction_by_epoch"],axis=0)
        r=np.sum(shape["q2_weighted_count_fraction_by_epoch"],axis=0)
        assert len(w)==len(self.z) and np.all(w>=0.) and np.all(r>=0.)
        assert np.all((w>0.)| (r==0.))
        self.w=w/w.sum()
        self.x=np.divide(r,w,out=np.zeros_like(r),where=w>0.)
        k=w[:-1]>0.
        assert np.all(self.x[:-1][k]>=self.z[:-1][k]*(1-1e-6))
        assert np.all(self.x[:-1][k]<=self.z[1:][k]*(1+1e-6))
        assert w[-1]==0. or self.x[-1]>=self.z[-1]*(1-1e-6)
        self.b=np.clip((self.x[:-1]-self.z[:-1])/np.diff(self.z),0.,1.)
        self.realized_second=np.sum(r)/np.sum(w)

    def recoil(self,u):
        K=self.clock
        ent=K.s(u)
        scale=TF**2*(ent/K.s(np.log(TF)))**(2/3)
        normal=self.omega/mp/self.fac*ent/(3*MPL**2*K.H(u)**2)
        p,z=self.x*scale,self.z*scale
        E=np.hypot(mp,np.sqrt(p))
        e=np.hypot(mp,np.sqrt(z))
        out=[]
        for m0,m1 in ((e,E),(z/(3*e),p/(3*E))):
            low=m0.copy()
            low[:-1]+=self.b*np.diff(m0)
            v=(normal*np.array([self.w@low,self.w@m1]))
            assert 0.<v[0]<=v[1]*(1.+2e-13)
            out.append(v)
        return np.array(out)

    def joint(self,u):
        r,p,_=self.prompt.ratio(u)
        ref=np.array([r.sum(axis=0),p.sum(axis=0)])
        before=np.array([r[1],p[1]])
        add=self.recoil(u)
        # Use the same stable spectrum in both comparisons. Carry over
        # the known correlation rather than subtracting unrelated bounds.
        delta=add-self.fraction*before[:,::-1]
        net=np.array([r[0]+(1-self.fraction)*r[1],p[0]+(1-self.fraction)*p[1]])+add
        H0=np.sqrt(1.+ref[0])
        H1=np.sqrt(1.+net[0])
        h=(delta[0,:,None]/(H0+H1)[None,:]).ravel()
        t=(-delta[0,:,None]/(H0*H1*(H0+H1))[None,:]).ravel()
        return ref,net,delta,add,np.array([h.min(),h.max()]),np.array([t.min(),t.max()])

    def discrete_control(self,T):
        # Independently realize a positive count/second-moment measure.
        # Put 37% at each bin mean and split the rest between finite edges
        # (use a finite high endpoint for the last open bin).  Compute
        # stress directly without calling either reduction inside joint().
        u=np.log(T)
        K=self.clock
        sr=TF**2*(K.s(u)/K.s(np.log(TF)))**(2/3)
        scale=K.s(u)/(3*MPL**2*K.H(u)**2)
        def integrate(mass,w,x,low,omega):
            high=np.r_[low[1:],max(16.*low[-1],32.*x[-1],1.)]
            ok=w>0.
            W,L,U,M=w[ok],low[ok],high[ok],x[ok]
            f=np.clip((M-L)/(U-L),0.,1.)
            q=np.r_[M,L,U]
            d=np.r_[.37*W,.63*W*(1.-f),.63*W*f]
            assert abs(d.sum()-1.)<1e-12
            assert abs((d@q)/(W@M)-1.)<1e-6
            E=np.hypot(mass,np.sqrt(sr*q))
            return scale*omega/(mass*self.fac)*np.array([d@E,d@(sr*q/(3*E))])
        F=self.prompt
        pieces=[integrate(a,F.Y[j],F.x[j],F.edges2,F.yield_count[j]*a*self.fac) for j,a in enumerate(F.mass)]
        extra=integrate(mp,self.w,self.x,self.z,self.omega)
        delta=extra-self.fraction*pieces[1]
        old=pieces[0]+pieces[1]
        new=old+delta
        h=delta[0]/(np.sqrt(1.+new[0])+np.sqrt(1.+old[0]))
        low,high=self.joint(u)[-2]
        assert low<=h<=high
        return dict(T_GeV=T,exact_positive_discrete_delta_H_over_H_radiation=float(h),
                    exact_positive_discrete_delta_stress_over_radiation=delta.tolist())

    def continuity(self,T):
        K=self.clock
        u=np.log(T)
        def terms(v):
            a,b,_,d,*_=self.joint(v)
            return np.array([a,b,d])*3*MPL**2*K.H(v)**2
        ref=terms(u)
        ans=[]
        for h in (2e-5,1e-5):
            derivative=(terms(u+h)[:,0]-terms(u-h)[:,0])/(2*h)
            expected=(3+K.grad(u))*ref.sum(axis=1)
            err=float(np.max(abs(derivative/expected-1.)))
            assert err<4e-7
            ans.append(err)
        return ans

    def integrated(self,T,n):
        K=self.clock
        u0,u1=np.log([T,TF])
        e=np.unique(np.clip(np.r_[u0,u1,K.bg.sg.x,np.linspace(u0,u1,50)],u0,u1))
        u,w=roots_legendre(n)
        a,b=e[:-1,None],e[1:,None]
        t=((a+b)+(b-a)*u)/2
        x=(b-a)*w/2
        out=np.zeros(3)
        for v,dv in zip(t.ravel(),x.ravel()):
            weight=dv*(1+K.grad(v)/3)/K.H(v)
            dt=self.joint(v)[-1]
            assert weight>0. and 0.<dt[0]<=dt[1]
            out+=weight*np.r_[1.,dt]
        return out


def run(eos,dist,shape,audit):
    D,S,A=[json.loads(p.read_text()) for p in (dist,shape,audit)]
    h=sha(eos)
    assert D["eos_sha256"]==S["eos_sha256"]==h and D["anchor_temperature_GeV"]==TF
    assert sha(dist) in A["record_sha256"].values() and S["reference_sha256"]==D["reference_sha256"]
    clock=SourceFold(eos,(4,9,8),65,4)
    fac=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    e,d=[np.asarray(x) for x in (D["boundary_p_over_T"],S["recoil_p_over_TF_lower_edges"])]
    sec=constants.hbar/(1e9*constants.e)
    results=[]
    for p,q,a in zip(D["cases"],S["episodes"],A["cases"],strict=True):
        assert p["T_R_GeV"]==q["T_R_GeV"]==a["test_root_GeV"] and p["Tmax_GeV"]==q["Tmax_GeV"]
        F=Shifted(clock,e,p,d,q["modes"][-1],a,fac)
        local=[]
        for T in (.001,.0007,.0001,.00003):
            old,new,delta,extra,dh,dt=F.joint(np.log(T))
            assert dh[1]<0. and dt[0]>0.
            local.append(dict(T_GeV=T,all_prompt_stress_over_radiation=old.tolist(),converted_stress_over_radiation=new.tolist(),
                              correlated_delta_stress_over_radiation=delta.tolist(),recoil_stress_over_radiation=extra.tolist(),
                              delta_H_over_H_radiation=dh.tolist(),delta_inverse_H_times_H_radiation=dt.tolist()))
        integ=[]
        for T in (.0007,.0001,.00003):
            low=F.integrated(T,5)
            hi=F.integrated(T,8)
            err=float(max(abs(low/hi-1.)))
            assert err<4e-7
            integ.append(dict(lower_T_GeV=T,fractional_elapsed_time_increase_over_radiation_elapsed=(hi[1:]/hi[0]).tolist(),
                              absolute_elapsed_time_increase_seconds=(sec*hi[1:]).tolist(),five_vs_eight_relative_defect=err))
        results.append(dict(ratio=a["ratio"],T_R_GeV=p["T_R_GeV"],removed_stable_count_fraction=F.fraction,
            scalar_recoil_second_moment_p_T=F.realized_second,probes=local,integrals=integ,
            continuity_relative_defects=[F.continuity(T) for T in (.00071,.000123,.000031)],
            discrete_controls=[F.discrete_control(T) for T in (.001,.0007,.0001,.00003)]))
    return dict(cases=results,record_hashes={p.name:sha(p) for p in (eos,dist,shape,audit)},stress_row_order=["energy","pressure"],
         scope="Correlated positive-bin bounds for homogeneous post-1 MeV expansion on the same radiation EOS, after removing an equal stable count on an unchanged stable shape and adding early spontaneous recoil. No post-anchor pair conversions, prior stress feedback, entropy change, measured nuclear yields or fitted Delta N_eff. These are not bounds on a full coupled model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","dist","shape","audit"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.dist,a.shape,a.audit),indent=2,allow_nan=False))
