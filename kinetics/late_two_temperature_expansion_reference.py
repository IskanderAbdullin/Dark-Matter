"""Independent FD/BE thermal moment closure with finite elastic heat."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from continuum_dark_normal_capacity_reference import critical_T,from_row,fugacity,moments,rule,solve
from continuum_dark_dilution_limits_reference import initial_shells
from late_signed_scalar_heat_reference import heat,nodes
from three_species_decay_onsager_reference import m,mp
from sm_dark_entropy_exchange_reference import StandardModel,MPL


def single(T,n,j):
    mass,bose=((mp,True),(m,False))[j]
    et=fugacity(T,mass,n,bose,20)
    _,K=moments(T,mass,et,bose)
    entropy=moments(T,mass,et,bose,entropy=True)
    return dict(T=T,eta=et,K=K,S=entropy)


def pressure(st,j):
    mass,bose=((mp,True),(m,False))[j]
    T,eta=st["T"],st["eta"]
    k,E,p,W=nodes(T,mass,18,64.)
    x=k/T-eta
    f=np.exp(-x)/(-np.expm1(-x)) if bose else expit(-x)
    return float(np.sum(W*p*p/E*f)*(1 if bose else 4)/3)


def inverse(n,target,j,quantity,guess):
    mass,bose=((mp,True),(m,False))[j]
    Tc=critical_T(n,20) if bose else 0.
    lo=max(guess/6.,Tc*(1+1e-13))
    hi=max(guess*5.,lo*2.)
    def func(u):
        return single(np.exp(u),n,j)[quantity]/target-1.
    while func(np.log(lo))>0.:
        lo*=.5
        if lo<=Tc: 
            lo=Tc*(1+1e-13)
            if func(np.log(lo))>0.: raise ValueError("specified thermal scalar state needs a non-normal branch")
        if lo<1e-24: raise ValueError("unrealizable closure moment")
    while func(np.log(hi))<0.: hi*=2.
    T=np.exp(brentq(func,np.log(lo),np.log(hi),xtol=3e-13))
    return single(T,n,j)


def common_entropy(n,S,guess):
    Tc=critical_T(n[0],20)
    def f(T):
        r=[single(T,x,j) for j,x in enumerate(n)]
        return sum(v["S"] for v in r),r
    lo=Tc*(1+1e-13)
    if f(lo)[0]>S: raise ValueError("common entropy branch not normal")
    hi=max(guess,lo*2.)
    while f(hi)[0]<S: hi*=2.
    T=np.exp(brentq(lambda q:f(np.exp(q))[0]/S-1.,np.log(lo),np.log(hi)))
    return f(T)[1]


def redshift(N,a,b,state):
    if a==b: return state,np.zeros(2)
    s=np.array([v["S"] for v in state])*a**3
    la,lb=np.log([a,b])
    out=[inverse(N[j]/b**3,s[j]/b**3,j,"S",state[j]["T"]*a/b) for j in range(2)]
    g,w=rule(5)
    trace=np.zeros(2)
    for x,y in zip((la+lb)/2+(lb-la)/2*g,(lb-la)/2*w):
        d=np.exp(x)
        for j in range(2):
            v=inverse(N[j]/d**3,s[j]/d**3,j,"S",state[j]["T"]*a/d)
            trace[j]-=3*y*d**3*pressure(v,j)
    delta=np.array([v["K"] for v in out])*b**3-np.array([v["K"] for v in state])*a**3
    assert max(abs(trace-delta)/(np.array([v["K"] for v in state])*a**3))<4e-8
    return out,trace


def kick(n,state,h,mode,orders):
    if mode=="off": return state,0
    total=state[0]["K"]+state[1]["K"]
    old=state[0]["K"]
    common=solve(n,total)["common_T_GeV"]
    ends=[single(common,n[j],j) for j in range(2)]
    shift=ends[0]["K"]-old
    if abs(shift)<1e-12*total: return state,0
    if mode=="equilibrate":
        return ends,0
    called=0
    def root(z):
        nonlocal called
        dst=old+z*shift
        s0=inverse(n[0],dst,0,"K",state[0]["T"])
        s1=inverse(n[1],total-dst,1,"K",state[1]["T"])
        q=heat(s0["T"],s1["T"],s0["eta"],s1["eta"],orders)
        if z==0.: assert q*shift>=0.
        called+=1
        return z-h*q/shift
    z=brentq(root,0.,1.,xtol=2e-9)
    return [inverse(n[0],old+z*shift,0,"K",state[0]["T"]),
            inverse(n[1],total-old-z*shift,1,"K",state[1]["T"])],called


def run(row,snap,eos,nx,start=3.,stop=410.,mode="full",orders=(5,5,6,3),initial="common"):
    cap=from_row(row,eos,20)
    mass_p=initial_shells(snap)
    N=np.array([np.sum(w) for _,_,w in mass_p])
    K=np.array([np.sum(w*p*p/(np.hypot(p,mass)+mass)) for mass,p,w in mass_p])
    assert abs(K.sum()/cap['retained_dark_kinetic_energy_GeV4']-1.)<5e-12
    if initial=="common":
        a=common_entropy(N/start**3,cap["entropy_density_GeV3"]/start**3,cap["common_T_GeV"]/start)
    else:
        a0=[inverse(N[j],K[j],j,"K",cap["common_T_GeV"]) for j in range(2)]
        a,_=redshift(N,1.,start,a0)
    Kin=np.array([v["K"] for v in a])*start**3
    S0=sum(v["S"] for v in a)*start**3
    pres=np.zeros(2); Q=np.zeros(2); dS=0.; ncall=0
    x=np.linspace(np.log(start),np.log(stop),nx+1)
    steps=[]
    tv=cap["actual_visible_T_GeV"]
    for l,r in zip(x,x[1:]):
        lo,mid,hi=np.exp([l,(l+r)/2,r])
        a,p=redshift(N,lo,mid,a);pres+=p
        H=np.sqrt(np.pi**2/45)*(tv/mid)**2/MPL
        b,nq=kick(N/mid**3,a,(r-l)/H,mode,orders)
        qq=np.array([w["K"]-v["K"] for w,v in zip(b,a)])*mid**3
        assert abs(qq.sum())<3e-10*Kin.sum()
        ds=sum(w["S"]-v["S"] for w,v in zip(b,a))*mid**3
        assert ds>-3e-11*S0
        dS+=ds;Q+=qq;ncall+=nq
        a,p=redshift(N,mid,hi,b);pres+=p
        n=N/hi**3
        cap_phi=moments(a[0]["T"],mp,0.,True)[0]
        steps.append(dict(b=float(hi),T_GeV=[v["T"] for v in a],scalar_density_over_saturation=float(n[0]/cap_phi)))
    Kf=np.array([v["K"] for v in a])*stop**3
    Eerror=max(abs(Kf-Kin-pres-Q))/Kin.sum()
    Serror=abs(sum(v["S"] for v in a)*stop**3-S0-dS)/S0
    assert max(Eerror,Serror)<4e-8
    return dict(nsteps=nx,mode=mode,initial=initial,start=start,stop=stop,
                quadrature=orders,heat_evaluations=ncall,
                kinetic_energy_ledger_defect=float(Eerror),entropy_ledger_defect=float(Serror),
                accumulated_scalar_comoving_heat_over_initial_total_comoving_kinetic_energy=float(Q[0]/Kin.sum()),
                fractional_entropy_production=float(dS/S0),steps=steps)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("run");p.add_argument("snapshot")
    p.add_argument("--steps",type=int,nargs="+",default=[12,24])
    p.add_argument("--start",type=float,default=3.)
    p.add_argument("--stop",type=float,default=410.)
    p.add_argument("--mode",choices=("off","equilibrate","full"),default="full")
    p.add_argument("--initial",choices=("common","isolated"),default="common")
    p.add_argument("--orders",type=int,nargs=4,default=(5,5,6,3))
    a=p.parse_args()
    row=json.loads(Path(a.run).read_text())["rows"][-1]
    snap=json.loads(Path(a.snapshot).read_text())["snapshots"][-1]
    assert row["nsteps"]==snap["steps"]
    eos=StandardModel(a.eos)
    print(json.dumps(dict(rows=[run(row,snap,eos,n,a.start,a.stop,a.mode,a.orders,a.initial) for n in a.steps]),indent=2,allow_nan=False))
