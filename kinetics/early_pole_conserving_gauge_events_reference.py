"""Energy-only gauge mean-field events, with reversible quantum occupation weights."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit



class Model:
    def __init__(self,g):
        self.g=g
        self.labels=['gamma','P','b','scalar']
        self.deg=np.array([2.,4.,4.,1.])
        self.p=[np.array([.35,.55,.7,.9,1.2,1.4,1.8,2.5]),
                np.array([.14,.34,.65,.86,1.11,1.35,1.44,1.68,2.]),
                np.array([.11,.33,.49,.64,.81,1.02,1.14,1.32,1.4389,1.445,1.56,1.75,2.2]),
                np.array([4.1e-6,.03,.15,.3])]
        self.M=np.array([0.,.05,.003,.00086724])
        E=[np.hypot(v,m) for v,m in zip(self.p,self.M)]
        self.E=np.concatenate(E)
        W=[]
        for j,p in enumerate(self.p):
            if j==3:
                W.append(np.array([.004,.006,.008,.01]))
            else:
                edges=np.r_[(p[0]*.6),(.5*(p[1:]+p[:-1])),p[-1]+.2]
                W.append(self.deg[j]*np.diff(edges**3)/(6*np.pi**2))
        self.slices={}
        b=0
        for j,k in enumerate(self.labels):
            self.slices[k]=slice(b,b+len(W[j]))
            b+=len(W[j])
        self.w=np.concatenate(W)
        self.stat=np.r_[np.ones(len(W[0])), -np.ones(len(W[1])+len(W[2])),np.ones(len(W[3]))]
        J=np.zeros((3,b))
        for j,k in enumerate(self.labels[:3]):
            sl=self.slices[k]
            J[j,sl]=2*np.pi**2/(self.deg[j]*self.E[sl])
        self.K=(g*g/np.pi**4)*J.T@np.array([[0,1,1],[1,1,0],[1,0,1]])@J
        self.charge=np.r_[np.ones(b-len(W[-1])),np.zeros(len(W[-1]))]
        self.topo=[('P',6,'scalar',0,'b',1.),
                   ('P',5,'scalar',1,'b',.5),
                   ('P',7,'scalar',2,'b',.65)]
        self.pair=[(2,4,'b',3,.22),(1,6,'P',4,.18)]
        f=np.concatenate([1/np.expm1(E[0]+.3),expit(-E[1]+.2),
                          expit(-E[2]/.55-1.4),1/np.expm1(E[3]+1.1)])
        self.init=self.w*f

    def ene(self,f):
        return self.E@f+.5*f@self.K@f

    def dE(self,f):
        return self.E+self.K@f

    def index(self,k,i):
        return self.slices[k].start+i

    def project(self,E,fix,site,kappa):
        vec=np.zeros_like(E)
        for i,v in fix:vec[i]+=v
        target=-E@vec
        sl=self.slices[site]
        j=np.searchsorted(E[sl],target)-1
        assert 0<=j<len(self.p[self.labels.index(site)])-1
        i=sl.start+j
        f=(E[i+1]-target)/(E[i+1]-E[i])
        assert 0<=f<=1
        vec[i]+=f;vec[i+1]+=1-f
        assert abs(E@vec)<1e-13
        assert abs(self.charge@vec)<1e-13
        return vec,kappa

    def events(self,E):
        out=[self.project(E,[(self.index(p,i),-1.),(self.index(q,j),1.)],b,k)
             for p,i,q,j,b,k in self.topo]
        for i,j,f,l,k in self.pair:
            out.append(self.project(E,[(self.index('gamma',i),-1.),(self.index('gamma',j),-1.),
                                       (self.index(f,l),1.)],f,k))
        return out

    def rates(self,n,balanced):
        f=n/self.w
        assert min(f)>0. and min(1+self.stat*f)>0.
        E=self.dE(n) if balanced else self.E
        logf=np.log(f)
        ln=np.log1p(self.stat*f)
        d=np.zeros_like(n)
        prod=0.
        for v,k in self.events(E):
            p=np.maximum(v,0.)
            m=np.maximum(-v,0.)
            x=m@logf+p@ln
            y=p@logf+m@ln
            gain=k*np.exp(y)*np.expm1(x-y)
            d+=v*gain
            prod+=gain*(x-y)
        assert prod>=0.
        return d,prod

    def equilibrium(self):
        mu=-.24*self.charge
        y=self.init.copy()
        for n in range(35):
            a=(self.dE(y)-mu)/.8
            nz=self.w/(np.exp(a)-self.stat)
            if max(abs(nz-y))<1e-13:
                y=nz
                break
            y=nz
        res=max(abs(self.rates(y,True)[0]))
        old=max(abs(self.rates(y,False)[0]))
        assert res<2e-12
        return dict(conservative_equilibrium_residual=res,bare_event_on_dressed_equilibrium_residual=old)

    def one_exact_finite_step(self):
        n=self.init; E=self.dE(n)
        v,_=self.events(E)[0]
        mask=np.flatnonzero(v>0)
        jr=[p for p in mask if self.slices['b'].start<=p<self.slices['b'].stop]
        assert len(jr)==2
        def z(f):
            u=v.copy()
            u[jr]=np.r_[f,1-f]
            return u
        dose=1e-4
        def defect(f):
            v=z(f)
            return self.dE(n+.5*dose*v)@v
        w=brentq(defect,0.,1.,xtol=1e-14)
        a=n+dose*z(w)
        b=n+dose*v
        return dict(exact_midpoint_error=self.ene(a)-self.ene(n),
                    uncorrected_tangent_error=self.ene(b)-self.ene(n),
                    weight_change=w-v[jr[0]])

    def simulate(self,balanced):
        ini=self.init.copy()
        dt,prod=self.rates(ini,balanced)
        def rhs(t,y):
            d,S=self.rates(y[:-1],balanced)
            return np.r_[d,S]
        sol=solve_ivp(rhs,(0.,.12),np.r_[ini,0.],method='DOP853',rtol=1e-11,atol=1e-14)
        assert sol.success
        end=sol.y[:-1,-1]
        def entropy(n):
            f=n/self.w
            b=1+self.stat*f
            return self.w@(self.stat*b*np.log(b)-f*np.log(f))
        r=dict(balanced=balanced,first_energy_derivative=self.dE(ini)@dt,
               energy_drift=self.ene(end)-self.ene(ini),gauge_count_drift=self.charge@(end-ini),
               entropy_integral=sol.y[-1,-1],entropy_change=entropy(end)-entropy(ini))
        assert abs(r['gauge_count_drift'])<1e-13
        assert abs(r['entropy_integral']-r['entropy_change'])<1e-10
        if balanced:
            assert abs(r['energy_drift'])<2e-13
        return r


def run():
    rows=[]
    for g in (0.,.005,.05):
        M=Model(g)
        rows.append(dict(g=g,event_step=M.one_exact_finite_step(),equilibrium=M.equilibrium(),data=[M.simulate(x) for x in (False,True)]))
    return dict(scope='Fixed-scale discrete energy-only reversible channel witness, not physical decay '
                      'geometry, photon kinetics, or FRW. Uses external fixed scalar mass.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
