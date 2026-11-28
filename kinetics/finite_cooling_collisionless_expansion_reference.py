"""Fixed-temperature homogeneous expansion bound for represented prompt PSDs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_portal_reference import TF, SourceFold
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m, mp


def moments(case,edges):
    Y=np.sum(case["count_fraction_per_epoch_species_bin"],axis=0)
    X=np.sum(case["p_over_T_squared_fractional_number_moment_per_epoch_species_bin"],axis=0)
    assert np.all(Y>=0.) and np.all(X>=0.)
    local=np.divide(X,Y,out=np.zeros_like(X),where=Y>0.)
    lo=edges**2
    ok=Y[:,:-1]>0
    assert np.all(local[:,:-1][ok]>=np.broadcast_to(lo[:-1],ok.shape)[ok]*(1-1e-6))
    assert np.all(local[:,:-1][ok]<=np.broadcast_to(lo[1:],ok.shape)[ok]*(1+1e-6))
    weight=np.clip((local[:,:-1]-lo[:-1][None,:])/np.diff(lo),0.,1.)
    return Y/Y.sum(axis=1)[:,None],local,weight


class Fluid:
    def __init__(self,K,edges,case,unit):
        self.K=K
        self.mass=np.array([mp,m])
        self.Y,self.x,self.W=moments(case,edges)
        self.edges2=edges**2
        self.yield_count=np.array([r["Omega_h2"] for r in case["species"]])/(unit*self.mass)

    def ratio(self,u):
        entropy=self.K.s(u)
        rel2=TF**2*(entropy/self.K.s(np.log(TF)))**(2/3)
        counts=self.yield_count*entropy
        m0=self.mass[:,None]
        rho_vis=3*MPL**2*self.K.H(u)**2
        zz=[]
        # Kinetic and pressure integrands both monotonic, concave in p².
        funcs=(lambda x,rho:x/(rho+m0),lambda x,rho:x/(3*rho))
        l=self.edges2[None,:]*rel2
        z=self.x*rel2
        for f in funcs:
            edges=f(l,np.sqrt(m0*m0+l))
            chord=edges.copy()
            chord[:,:-1]+=self.W*np.diff(edges,axis=1)
            upper=f(z,np.sqrt(m0*m0+z))
            fr=np.array([(self.Y*chord).sum(axis=1),(self.Y*upper).sum(axis=1)]).T*counts[:,None]/rho_vis
            assert np.all(fr[:,0]<=fr[:,1]*(1+1e-7))
            zz.append(fr)
        rest=counts*self.mass/rho_vis
        rho=zz[0]+rest[:,None]
        p=zz[1]
        return rho,p,rest

    def continuity_check(self,T):
        u=np.log(T)
        def rho(v):
            d=self.ratio(v)[0][:,1]
            return d*3*MPL**2*self.K.H(v)**2
        p=self.ratio(u)[1][:,1]*3*MPL**2*self.K.H(u)**2
        derivative=(rho(u+1e-5)-rho(u-1e-5))/(2e-5)
        predicted=(3+self.K.grad(u))*(rho(u)+p)
        err=np.max(abs(derivative/predicted-1.))
        assert err<4e-7
        return float(err)

    def clock(self,T,n):
        lo,hi=np.log([T,TF])
        x,w=roots_legendre(n)
        e=np.unique(np.clip(np.r_[lo,hi,self.K.bg.sg.x,np.linspace(lo,hi,50)],lo,hi))
        u=((e[:-1,None]+e[1:,None])+np.diff(e)[:,None]*x)/2
        v=np.diff(e)[:,None]*w/2
        a=np.zeros(3)
        for uu,b in zip(u.ravel(),v.ravel()):
            t=b*(1+self.K.grad(uu)/3)/self.K.H(uu)
            d=self.ratio(uu)[0].sum(axis=0)
            f=np.sqrt(1+d)
            a+=t*np.r_[1.,d/(f*(1+f))]
        return a


def run(eos,dist):
    reference=json.loads(Path(dist).read_text())
    esh=hashlib.sha256(Path(eos).read_bytes()).hexdigest()
    assert esh==reference["eos_sha256"]
    bg=SourceFold(eos,(4,9,8),65,4)
    conv=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    assert TF==reference["anchor_temperature_GeV"]
    edge=np.asarray(reference["boundary_p_over_T"])
    ts=constants.hbar/(1e9*constants.e)
    R=[]
    for case in reference["cases"]:
        F=Fluid(bg,edge,case,conv)
        rows=[]
        for T in (.001,.0007,.0001,.00003,.00001):
            d,p,r=F.ratio(np.log(T))
            n=d.sum(axis=0)
            delta=n/(np.sqrt(1+n)+1.)
            rows.append(dict(T_GeV=T,rho_dark_species_over_Rho_SM=d.tolist(),
                P_dark_species_over_Rho_SM=p.tolist(),rho_dark_rest_over_Rho_SM=r.tolist(),
                fractional_H_increase=delta.tolist()))
        times=[]
        for T in (.0007,.0001,.00003):
            a=F.clock(T,8);b=F.clock(T,5)
            err=np.max(abs(a/b-1))
            assert err<1e-7
            times.append(dict(T_end_GeV=T,visible_only_elapsed_seconds=float(a[0]*ts),
                fractional_elapsed_time_shortening=(a[1:]/a[0]).tolist(),max_5_vs_8_relative_defect=float(err)))
        R.append(dict(T_R_GeV=case["T_R_GeV"],Tmax_GeV=case["Tmax_GeV"],
            fixed_T=rows,time_integrals=times,continuity_check_at_temperatures_GeV=[
                dict(T_GeV=T,relative_energy_pressure_continuity_defect=F.continuity_check(T))
                for T in (.00071,.000123,.000031)]))
    return dict(eos_sha256=esh,distribution_sha256=hashlib.sha256(Path(dist).read_bytes()).hexdigest(),
        scope="Post-1MeV collisionless stress and homogeneous clock on the specified visible-radiation background with and without the positive represented dark energy. No standard-matter double count in the radiation-only comparator; no nuclear abundances, coupled history, dark production before the anchor, or DeltaNeff inference.",cases=R)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("distribution")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.distribution),indent=2,allow_nan=False))
