"""Local FD dark-temperature drift, distinct from net heat or pair production."""
import json
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize_scalar, root

from boost_axis_elastic_heat_reference import integrate, leptons
from fd_number_heat_state_reference import evaluate, state
from positive_fd_portal_quadrature_reference import PositivePortal


class Control:
    def __init__(self, fine):
        self.portal=PositivePortal((13,36,26) if fine else (9,24,16),tail=65. if fine else 55.)
        self.order=22 if fine else 16
        self.cut=35. if fine else 29.
        self.nodes=174 if fine else 112

    @lru_cache(maxsize=2048)
    def flow(self,T,r,eta,elastic=True):
        D=T*r
        F,Qt,S=self.portal.flow(T,D,eta)
        Q=Qt-sum(integrate(T,D,eta,m,n=self.order,cut=self.cut) for m in leptons) if elastic else Qt
        x=evaluate(D,eta,nodes=self.nodes)
        chi=x['chemical_energy']
        G=Q-2*chi*F
        sigma=2*eta*F+(D-T)/(D*T)*Q
        if sigma<0. or x['cv_fixed_number']<=0.:
            raise ValueError("Entropy or fixed-density heat-capacity sign.")
        return dict(G=float(G),F=float(F), Q=float(Q), Qt=float(Qt),
                    sigma=float(sigma),chi=float(chi), n=float(x['n']),
                    rho=float(x['rho']),cv=float(x['cv_fixed_number']),timelike_entropy=float(S))

    def scan(self,T,r,maximum=5.,cells=100,elastic=True):
        f=lambda eta: self.flow(T,r,eta,elastic)['G']
        grid=np.linspace(-maximum,maximum,cells+1)
        vals=np.array([f(float(a)) for a in grid])
        brackets=np.flatnonzero(vals[:-1]*vals[1:]<0.)
        return [brentq(f,float(grid[i]),float(grid[i+1]),xtol=1e-12) for i in brackets]

    def fold(self,T,elastic=True):
        # Track the local positive-eta maximum in a predeclared rectangle.
        def peak(r):
            f=lambda eta:self.flow(T,r,float(eta),elastic)['G']
            m=minimize_scalar(lambda z:-f(z),bounds=(1.4,4.6),method="bounded",
                              options=dict(xatol=2e-8))
            if not m.success or not 1.5<m.x<4.5:
                raise ValueError("Peak was not interior to the local rectangle.")
            return float(m.x),f(float(m.x))
        bounds={(.015,True):(.84,.86),(.015,False):(.70,.75),(.1,True):(.80,.82),(.1,False):(.63,.67)}
        r=brentq(lambda z:peak(z)[1],*bounds[T,elastic],xtol=8e-12)
        eta,G=peak(r)
        base=self.flow(T,r,eta,elastic)
        norm=max(abs(base['Q']),abs(2*base['chi']*base['F']))
        f=lambda a,b:self.flow(T,a,b,elastic)['G']/norm
        h=1e-4
        ratio_slope=(f(r+h,eta)-f(r-h,eta))/(2*h)
        curvature=(f(r,eta+h)+f(r,eta-h)-2*f(r,eta))/(h*h)
        eta_grad=(f(r,eta+h)-f(r,eta-h))/(2*h)
        if not (abs(G/norm)<2e-9 and abs(eta_grad)<3e-6 and ratio_slope>0 and curvature<0):
            raise ValueError(f"Not a local max fold: {T,r,eta,G/norm,eta_grad,ratio_slope,curvature}.")
        rr=root(lambda v:[f(*v),(f(v[0],v[1]+h/5)-f(v[0],v[1]-h/5))/(2*h/5)],
                [r+3e-4,eta-.01],tol=5e-9)
        if not rr.success or np.max(abs(rr.x-[r,eta]))>1e-6:
            raise ValueError("Independent coupled root did not match fold.")
        ra,rb=r-5e-4,r+5e-4
        p1,p2=peak(ra),peak(rb)
        if not p1[1]<0<p2[1] or max(f(rb,1.4),f(rb,4.6))>=0.:
            raise ValueError("Fold sides did not bracket the local zero pair.")
        pair=[brentq(lambda z:f(rb,z),1.4,p2[0]),brentq(lambda z:f(rb,z),p2[0],4.6)]
        return dict(T=T,ratio=r,eta=eta,normalized_thermal_residual=G/norm,
                    ratio_slope=ratio_slope,eta_curvature=curvature,
                    normalized_peak_below=p1[1]/norm,normalized_peak_above=p2[1]/norm,
                    paired_zero_above=pair,independent_coupled_root=rr.x.tolist())

    def inversion(self,T,r,eta):
        D=T*r
        f=self.flow(T,r,eta)
        F,Q,G,cv,n,rho=[f[k] for k in ('F','Q','G','cv','n','rho')]
        natural=min(n/max(2*abs(F),1e-230),rho/max(abs(Q),1e-230))
        predicted=-G*natural/cv
        h=3e-5
        def temp(u):
            nr=n-2*u*natural*F
            er=rho-u*natural*Q
            return brentq(lambda x:state(x,nr,nodes=self.nodes)['rho']-er,.97*D,1.03*D,xtol=1e-15)
        a,b,c,d=[temp(j*h) for j in (-2,-1,1,2)]
        diff=(a-8*b+8*c-d)/(12*h)
        if abs(diff-predicted)>5e-8*D:
            raise ValueError("Finite number/energy-ledger temperature inversion failed.")
        return dict(eta=eta,analytic_temperature_drift=predicted,finite_ledger_temperature_drift=diff,**f)


def run():
    fine,coarse=Control(True),Control(False)
    folds=[]
    for T in (.015,.1):
        a,b=fine.fold(T),coarse.fold(T)
        if abs(a['ratio']-b['ratio'])>2e-5 or abs(a['eta']-b['eta'])>2e-4:
            raise ValueError("Fold moved under quadrature refinement.")
        folds.append(dict(fine=a,coarse=b,timelike_only=fine.fold(T,False)))
    cases=[]
    expected_counts=((.015,.8,0),(.015,.95,2),(.05,.8,2),(.1,.8,0),(.1,.95,2),
                     (.3,.8,2),(.6882487043,.8,1))
    for T,r,N in expected_counts:
        a,b=fine.scan(T,r,5.,120),coarse.scan(T,r,5.,64)
        if len(a)!=N or len(b)!=N or max([abs(x-y) for x,y in zip(a,b,strict=True)]+[0.])>2e-4:
            raise ValueError(f"Unexpected zeros in the reported eta interval: {T,r,a,b}.")
        cases.append(dict(T=T,ratio=r,eta_interval=[-5.,5.],fine_roots=a,coarse_roots=b,
                          response_checks=[fine.inversion(T,r,e+d) for e in a for d in (-.015,0.,.015)]))
    common=[]
    for T in (.015,.1,.3,.6882487043):
        _,_,_,M0,M1,M2,_=fine.portal.moments(T,T,0.)
        K=sum(integrate(T,T,0.,m,n=fine.order,cut=fine.cut,linear=True) for m in leptons)
        chi=evaluate(T,0.,nodes=fine.nodes)['chemical_energy']
        ex=-(M2/T+K-2*chi*M1/T)/(2*(M1-2*chi*M0))
        h=3e-4
        rr=[brentq(lambda z:fine.flow(T,1.+j*h,z)['G'],-.03,.03) for j in (-1,1)]
        fd=(rr[1]-rr[0])/(2*h)
        if abs(fd/ex-1)>8e-6:
            raise ValueError("Common temperature-root slope did not replay.")
        common.append(dict(T=T,analytic_eta_slope=ex,finite_eta_slope=fd))
    return dict(scope="Dark-temperature collision source at fixed instantaneous states; no expansion or trajectory.",
                folds=folds,finite_eta_cases=cases,common_ratio_slope=common)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
