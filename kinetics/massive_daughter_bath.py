"""Finite scalar mass in energy-closed 1<->2 kinetics, massless fermion daughter."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from energy_closed_bath import EnergyClosedBath


def allowed_area(elo,ehi,rlo,rhi,x,m):
    K=(x*x-m*m)/2
    assert K>0
    def positive_integral(E):
        p=np.sqrt(np.maximum((E-m)*(E+m),0.))
        low=np.maximum(rlo,K/(E+p))
        high=rhi if m==0 else np.minimum(rhi,K*(E+p)/(m*m))
        dr=np.maximum(high-low,0.)
        return E*dr-K/2*np.log1p(dr/low)-m*m/(4*K)*dr*(2*low+dr)
    return np.maximum(positive_integral(ehi)-positive_integral(elo),0.)


class MassiveDaughterBath(EnergyClosedBath):
    def __init__(self,n=200,mass_ratio=.001,**kwargs):
        super().__init__(n=n,**kwargs)
        assert 0<=mass_ratio<1
        self.mass_ratio=mass_ratio
        l,h=self.edges[:-1],self.edges[1:];u,w=roots_legendre(8)
        self.ps=(h[:,None]+l[:,None])/2+(h-l)[:,None]/2*u
        self.ws=(h-l)[:,None]/2*w*self.ps**2
        self.ws/=self.ws.sum(axis=1)[:,None]

    def geometry(self,x):
        if x==self.geom_x:return
        m=x*self.mass_ratio
        if m==0:
            super().geometry(x);return
        e=np.hypot(self.edges,m);energy=np.hypot(self.ps,m)
        self.energy_phi=np.sum(self.ws*energy,axis=1)
        tr=np.sum(self.ws*m*m/energy,axis=1)
        self.ew=np.r_[self.vol*self.energy_phi,4*self.vol*self.mean]
        self.trace_weights=np.r_[self.vol*tr,np.zeros(self.n)]
        self.pair_E=self.energy_phi[:,None]+self.mean[None,:]
        area=allowed_area(e[:-1,None],e[1:,None],self.lr,self.hr,x,m)
        # Gamma is the actual massive-channel rest width, eta=2*p_star/M.
        eta=1-self.mass_ratio**2
        self.R0=4*x/(2*np.pi**2)*area/eta*(self.kappa*x*x)
        self.geom_x=x

    def test_geometry(self):
        rng=np.random.default_rng(53)
        for x in [.03,.1,1.,5.]:
            m=x*self.mass_ratio;K=(x*x-m*m)/2;self.geometry(x)
            e=np.hypot(self.edges,m)
            max_err=0.
            # Sample cells around the lower-energy kinematic boundary.
            for i in rng.integers(0,self.n,80):
                j0=np.searchsorted(self.edges,K/(e[i]+self.edges[i]))-1
                for j in [j0-1,j0,j0+1,int(rng.integers(0,self.n))]:
                    if j<0 or j>=self.n:continue
                    a,b=self.edges[j:j+2];lo,hi=e[i:i+2]
                    def integrand(r):return max(hi-max(lo,(K/r+m*m*r/K)/2),0.)
                    pts=[]
                    for E in [lo,hi]:
                        pp=np.sqrt(max((E-m)*(E+m),0.));pts.append(K/(E+pp))
                        if m>0:pts.append(K*(E+pp)/(m*m))
                    direct=quad(integrand,a,b,points=[v for v in pts if a<v<b],epsabs=1e-18)[0]
                    computed=float(allowed_area(lo,hi,a,b,x,m))
                    err=abs(computed-direct)/max((b-a)*(hi-lo),1e-16)
                    max_err=max(max_err,err)
            assert max_err<2e-4,(x,max_err)
            # Exact equilibrium for discretized energy representatives.
            th=self.thermo(x,0.,0.);N=th['Neq'];y=np.r_[1/np.expm1(self.energy_phi),1/(np.exp(self.mean)+1),N,0.,0.]
            self.set_state(x,y)
            assert np.sum(abs(self.collision_rhs(y[:2*self.n])[:self.n])*self.vol)<1e-10*max(self.A.sum(),1e-30)


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=240);p.add_argument('--x0',type=float,default=.1)
    p.add_argument('--xend',type=float,default=40.);p.add_argument('--kappa',type=float,default=.1383)
    p.add_argument('--relaxation',type=float,default=100.);p.add_argument('--mass-ratios',type=float,nargs='+',default=[0.,1e-5,.001,.01,.03])
    p.add_argument('--rtol',type=float,default=2e-6);p.add_argument('--output',required=True);p.add_argument('--keep-psd',action='store_true');a=p.parse_args();rows=[]
    for mu in a.mass_ratios:
        m=MassiveDaughterBath(a.n,mass_ratio=mu,kappa=a.kappa,relaxation=a.relaxation)
        m.test_geometry()
        x=.2;z=np.log(x);y=np.r_[np.full(m.n,.1),np.full(m.n,.05),m.thermo(x,0.,0.)['Neq']*.8,.01,0.]
        rng=np.random.default_rng(17);d=rng.normal(size=len(y));d/=np.linalg.norm(d)
        eps=1e-6;j=m.jac_log(z,y)@d;fd=(m.rhs_log(z,y+eps*d)-m.rhs_log(z,y-eps*d))/(2*eps)
        assert np.linalg.norm(j-fd)/max(np.linalg.norm(j),1e-30)<1e-5
        r=m.evolve_closed(a.x0,a.xend,a.rtol,a.keep_psd);r['mass_ratio']=mu;rows.append(r)
        print(json.dumps({'mass_ratio':mu,'n':a.n,'nfev':r['nfev'],**r['final']}),flush=True)
    Path(a.output).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
