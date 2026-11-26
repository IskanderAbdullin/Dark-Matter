"""Higher-reheat scalar count threshold in the prescribed vacuum lepton source."""
import argparse
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion,proxy
from empty_lepton_pair_spectrum_reference import Source
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp


class Fold:
    def __init__(self,table,order,time_order,tf=.001,upper=.65):
        self.bg=TabulatedRadiation(table,.05)
        self.src=[Source(order,dark_mass=mm) for mm in (m,M)]
        self.x,self.w=roots_legendre(time_order)
        self.tf=tf
        a,b=np.log([tf,upper])
        self.edges=np.unique(np.r_[a,np.clip(self.bg.sg.x,a,b),b])
        self.sf,self.rf=self.densities(a)
        cumul=[np.zeros(4)]
        for l,r in zip(self.edges[:-1],self.edges[1:]):
            cumul.append(cumul[-1]+self.panel(l,r))
        self.cum=np.asarray(cumul)

    def densities(self,u):
        s=2*np.pi**2/45*np.exp(self.bg.entropy_log(u))
        rho=np.pi**2/30*np.exp(4*u+self.bg.rg(u))
        return s,rho

    def panel(self,lo,hi):
        out=np.zeros(4)
        for u,x in zip((lo+hi)/2+(hi-lo)/2*self.x,(hi-lo)/2*self.w):
            T=np.exp(u)
            s,rho=self.densities(u)
            H=np.sqrt(rho/3)/MPL
            dt=x*(3+self.bg.sg.derivative()(u))/(3*H)
            a=self.sf/s
            for j,src in enumerate(self.src):
                q,_,_,N,_=src.event_moments(T)
                out[j]+=2*N.sum()*dt/s
                E=np.sum(q*N)
                # Radiation-loss floor if energy stays in this isolated dark sector.
                out[2+j]+=E*dt*a**(4./3)/self.rf
        return out

    def sample(self,tr):
        u=np.log(tr)
        j=np.searchsorted(self.edges,u,side="right")-1
        assert 0<=j<len(self.edges)
        return self.cum[j]+self.panel(self.edges[j],u)


def report(table):
    f=conversion(table)["Omega_h2_per_GeV_times_yield"]
    trials=[]
    for order,ot in (((7,18,14),5),((10,27,20),7)):
        engine=Fold(table,order,ot)
        def sample(T):
            d,h,erd,erh=engine.sample(T)
            phi,chi=proxy([h,d+h],f)
            y=phi+chi*mp/m
            rho_dark_min=erd+erh
            # Ratio to radiation at the endpoint if stable particles survive.
            stable_rest=(d+h)*m*engine.sf/engine.rf
            return dict(T_reheat_GeV=T,scalar_upper_Omega_h2_from_total_count=y,
                         scalar_and_stable_fixed_source_Omega_h2=[phi,chi],
                         minimal_Redshifted_dark_energy_if_retained_over_visible_at_1MeV=rho_dark_min,
                         stable_rest_energy_over_visible_radiation_at_1MeV_if_count_preserved=stable_rest)
        rows=[sample(T) for T in (.05,.075,.10,.13,.16,.20,.25,.30,.40,.50,.60)]
        assert np.all(np.diff([q["scalar_upper_Omega_h2_from_total_count"] for q in rows])>0.)
        roots=[]
        for portion in (.5,.99,1.):
            target=.120*portion
            bracket=next((a,b) for a,b in zip(rows,rows[1:]) if
                         a["scalar_upper_Omega_h2_from_total_count"]<target<=b["scalar_upper_Omega_h2_from_total_count"])
            T=brentq(lambda T:sample(T)["scalar_upper_Omega_h2_from_total_count"]-target,
                bracket[0]["T_reheat_GeV"],bracket[1]["T_reheat_GeV"],xtol=1e-8)
            roots.append(dict(fraction_of_reference=portion,**sample(T)))
        trials.append(dict(order=order,time_order=ot,roots=roots,samples=rows))
    assert max(abs(a["T_reheat_GeV"]/b["T_reheat_GeV"]-1) for a,b in zip(trials[0]["roots"],trials[1]["roots"]))<5e-5
    return trials


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("table")
    a=p.parse_args()
    print(json.dumps(dict(trials=report(a.table),
         scope="Lepton-only, vacuum-final direct and prompt production on prescribed radiation. All-convert-to-scalar count ceiling; no annihilation, feedback or other channels used in the source. Energy marked as a floor only under pure retained dark injection without later visible losses; not BBN or a physical high-reheat solution."),
         indent=2,allow_nan=False))
