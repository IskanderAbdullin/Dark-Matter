"""Two-temperature trial with actual SM EOS and conserved massive partner number."""
import argparse
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

C=np.pi**2/30
MPL=2.435e18 # reduced Planck mass in GeV, H output only
MCHI=.02
Z,W=roots_legendre(152)
Y=(Z+1)*28
W=W*28*Y*Y*4/(2*np.pi**2)


class StandardModel:
    def __init__(self,file):
        t,gs,gr=np.loadtxt(file).T
        ok=(t>=.009)&(t<=1.1)
        u=np.log(t[ok])
        assert u[0]<np.log(.01) and u[-1]>0.
        self.gs=PchipInterpolator(u,gs[ok],extrapolate=False)
        self.gr=PchipInterpolator(u,gr[ok],extrapolate=False)
        self.slope=self.gs.derivative()
        t0=np.exp(u[0])
        p0=C*t0**4*(4*gs[ok][0]/3-gr[ok][0])
        self.press=solve_ivp(lambda x,y:[4*C/3*self.gs(x)*np.exp(4*x)],
                        (u[0],u[-1]),[p0],rtol=3e-11,atol=1e-19,
                        max_step=.006,dense_output=True)
        assert self.press.success
        self.lo,self.hi=u[[0,-1]]
    def state(self,T):
        x=np.log(T)
        assert self.lo<=x<=self.hi
        s=4*C/3*T**3*self.gs(x)
        P=self.press.sol(x)[0]
        rho=T*s-P
        cv=4*C/3*T**3*(3*self.gs(x)+self.slope(x))
        assert rho>0 and cv>0
        return float(rho),float(P),float(s),float(cv)


def dark(T,n):
    # Relativistic one-dof scalar; the massive fermion's number is fixed.
    E=T*np.hypot(Y,MCHI/T)
    w=T**3*W
    rho=C*T**4
    pressure=rho/3
    cv=4*rho/T
    s=4*rho/(3*T)
    if n==0.:
        return rho,pressure,s,cv,0.,None
    kappa=Y*Y/(np.hypot(Y,MCHI/T)+MCHI/T)
    f=lambda u:expit(u-kappa)
    rel_eta=brentq(lambda a:(w@f(a))/n-1.,-100.,150.,xtol=1e-12)
    eta=rel_eta+MCHI/T
    fn=f(rel_eta)
    assert abs((w@fn)/n-1)<3e-12
    r=w@(E*fn)
    p=w@((T*Y)**2/(3*E)*fn)
    wt=w*fn*(1-fn)
    dn=wt.sum()
    rn=(wt@E)/dn
    cvn=wt@((E-rn)**2)/(T*T)
    return rho+r,pressure+p,s+(r+p-eta*T*n)/T,cv+cvn,float(rn),float(eta)


def run_case(eos, coupling, td0_ratio=1.,heavy=True,rtol=2e-10,form="temperature",law=None):
    T0=.6882487043
    TD0=td0_ratio*T0
    n0=.015*TD0**3 if heavy else 0.
    r0,p0,ss0,_=eos.state(T0)
    d0,q0,sd0,*_=dark(TD0,n0)
    ss_end=eos.state(.01)[2]
    z0=np.r_[[ss0,sd0] if form=="entropy" else [np.log(T0),np.log(TD0)],0.,0.,0.,0.]
    def end(l,y):
        return (y[0]*np.exp(-3*l)-ss_end) if form=="entropy" else (y[0]-np.log(.01))
    def temperatures(l,y):
        if form=="temperature":
            return np.exp(y[:2])
        s,sd=y[:2]*np.exp(-3*l)
        Ts=brentq(lambda x:eos.state(np.exp(x))[2]/s-1.,eos.lo,eos.hi,xtol=1e-12)
        n=n0*np.exp(-3*l)
        if n==0.:
            Td=(sd/(4*C/3))**(1/3)
        else:
            guess=TD0*np.exp(-l)
            Td=brentq(lambda v:dark(v,n)[2]/sd-1.,.2*guess,5*guess,xtol=1e-14)
        return np.exp(Ts),Td
    end.terminal=True
    end.direction=-1
    def rhs(l,y):
        Ts,Td=temperatures(l,y)
        n=n0*np.exp(-3*l)
        rs,ps,ss,cs=eos.state(Ts)
        rd,pd,sd,cd,rn,eta=dark(Td,n)
        H=np.sqrt((rs+rd)/3)/MPL
        q=(coupling*(cs*cd/(cs+cd))*(Ts-Td) if law is None
           else law(Ts,Td,n)/H)
        ei=np.exp(4*l)
        prod=np.exp(3*l)*q*(1/Td-1/Ts)
        assert prod>=-1e-25
        pair=(-np.exp(3*l)*q/Ts,np.exp(3*l)*q/Td) if form=="entropy" else (
                     (-3*(rs+ps)-q)/(cs*Ts),(-3*(rd+pd-n*rn)+q)/(cd*Td))
        return np.r_[pair,
                     ei*(rs-3*ps), ei*(rd-3*pd),ei*q,prod]
    sol=solve_ivp(rhs,(0.,10.),z0,rtol=rtol,atol=rtol*2e-3,max_step=.025,
                 events=end,method="DOP853")
    assert sol.success and len(sol.t_events[0])==1
    le=sol.t[-1]
    y=sol.y[:,-1]
    Ts,Td=temperatures(le,y)
    ns=n0*np.exp(-3*le)
    rs,ps,ss,_=eos.state(Ts)
    rd,pd,sd,_,_,eta=dark(Td,ns)
    scale=np.exp(4*le)
    errors=[(scale*rs-r0)-(y[2]-y[4]),(scale*rd-d0)-(y[3]+y[4])]
    rel=max(abs(np.array(errors)))/(abs(r0)+abs(d0))
    sf=np.exp(3*le)*(ss+sd)
    se=abs(sf-ss0-sd0-y[5])/(ss0+sd0)
    assert rel<3e-6 and se<3e-6 and y[5]>=0.
    tdec=(float(eos.gs(np.log(Ts))/eos.gs(np.log(T0))))**(1/3)*td0_ratio
    if coupling==0. and not heavy and law is None:
        assert abs(Td/Ts/tdec-1.)<2e-7
    return dict(lambda_over_hubble=coupling if law is None else None,heavy=heavy,initial_dark_to_sm_T=td0_ratio,
                exact_initial_SM_T_GeV=T0,terminal_sm_T_GeV=Ts,terminal_dark_T_GeV=Td,
                terminal_temperature_ratio=Td/Ts,adiabatic_relativistic_ratio=tdec,
                terminal_log_scale_factor=le,terminal_fermion_number=ns,
                terminal_fermion_mu_over_T=eta,terminal_H_GeV=float(np.sqrt((rs+rd)/3)/MPL),
                comoving_heat_sm_to_dark=float(y[4]),entropy_produced=float(y[5]),
                relative_energy_ledger_error=float(rel),relative_entropy_ledger_error=float(se),
                relative_afferent_sm_entropy_change=float(np.exp(3*le)*ss/ss0-1),
                relative_afferent_dark_entropy_change=float(np.exp(3*le)*sd/sd0-1),
                nfev=sol.nfev,evolved_variables=form)


def run(table):
    eos=StandardModel(table)
    jobs=[(0.,1.,False),(0.,1.,True),(.3,1.,True),(3.,1.,True),(3.,1.15,True),(3.,.7,True)]
    rows=[run_case(eos,*a) for a in jobs]
    convergence=run_case(eos,3.,1.,True,2e-12)
    controls=[run_case(eos,*args,rtol=1e-9,form="entropy") for args in jobs[:4]]
    assert max(abs(x["terminal_temperature_ratio"]-y["terminal_temperature_ratio"]) for x,y in zip(rows,controls))<2e-6
    assert abs(convergence["terminal_temperature_ratio"]-rows[3]["terminal_temperature_ratio"])<2e-6
    return dict(rows=rows,tighter_control=convergence,entropy_variable_controls=controls,
                model="Phenomenological antisymmetric heat flow; not a calibrated SM channel or cosmological relic.")
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("table")
    print(json.dumps(run(ap.parse_args().table),indent=2,allow_nan=False))
