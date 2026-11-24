"""Energy moment upper on scalar multiplication with selected dark channels."""
import math,json
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import kve
from kinetics.phi_bose_portal_low_reheat_gate import MPL,M,GAMMA,GSTAR,S0,RHOC,MPHI,M1,M2,COUP2,LEPTONS,gamma
from kinetics.phi_bose_portal_vacuum_depletion_envelope import Icoef,make_table
YUK=.001
KINV=YUK**4/(math.pi*M2*M2)
BETA_MIN=math.sqrt(1-(M2/M1)**2)
B=(M1+M2)**2-MPHI**2
DELTA=M1**2-M2**2
KD=YUK**4/(16*math.pi*M1**2*BETA_MIN)*(1+B**2/DELTA**2)
JCOEF=288*COUP2/(math.pi**3*M**4)*MPL/(1.66*math.sqrt(GSTAR)*3)

def energy_source(T,mchi,ml=0.,cut=100.):
    u0=2*max(mchi,ml)/T
    def f(x):
        u=u0+x;s=(T*u)**2;a=ml*ml/s;b=mchi*mchi/s
        pref=(1+2*a)*(1+2*b)*math.sqrt(max(0,1-4*a))*math.sqrt(max(0,1-4*b))
        return u**6*pref*(2*kve(1,u)+u*kve(0,u))*math.exp(-x)/((s-M*M)**2+(M*GAMMA)**2)
    return COUP2*T**9/(48*math.pi**5)*math.exp(-u0)*quad(f,0,cut,epsabs=1e-30,epsrel=1e-8)[0]

def make_energy_table():
    ts=np.geomspace(.0005,.1,160);vals=[]
    for T in ts:
        s=2*math.pi**2/45*GSTAR*T**3;h=1.66*math.sqrt(GSTAR)*T*T/MPL
        vals.append(sum(energy_source(T,m,ml) for m in (M1,M2) for ml in LEPTONS)/(s*h))
    return PchipInterpolator(np.log(ts),np.log(vals))

def direct_sigma_v(s):
    if s<=4*M1*M1:return 0.
    bi=math.sqrt(1-4*M2*M2/s);bf=math.sqrt(1-4*M1*M1/s)
    A0=s/2-M2**2-M1**2;B0=s*bi*bf/2
    tp=-DELTA**2/(A0+B0);tm=-A0-B0
    zlo=MPHI**2-tp;zhi=MPHI**2-tm
    integral=(zhi-zlo)+2*B*math.log(zhi/zlo)+B*B*(1/zlo-1/zhi)
    return YUK**4*integral/(8*math.pi*s*s*bi)

def audit():
    check=[]
    for T in (.005,.02,.1):
        for m in (M1,M2):
            h=T*1e-5;d=T*T*(gamma(T+h,m)-gamma(T-h,m))/(2*h);v=energy_source(T,m)
            check.append(abs(v/d-1))
    tab=make_energy_table();ntab=make_table();rows=[]
    for TR in (.005,.008,.01,.015,.02,.025,.03,.05,.1):
        def src(z):
            T=TR*math.exp(-z)
            return math.exp(float(tab(math.log(T)))) if T>=.0005 else 0.
        X=quad(src,0,5,epsabs=1e-35,epsrel=1e-8)[0]
        def rhs(z,y):return [src(z)/X,math.exp(-z)*y[0]**2]
        sol=solve_ivp(rhs,(0,5),[0.,0.],rtol=2e-9,atol=1e-11)
        assert sol.success
        f=sol.y[1,-1]+math.exp(-5)*sol.y[0,-1]**2
        assert abs(sol.y[0,-1]-1)<1e-6 and 0<f<1
        Y_in=quad(lambda z:math.exp(float(ntab(math.log(TR)-z))) if TR*math.exp(-z)>=.0005 else 0.,0,5,epsabs=1e-35,epsrel=1e-7)[0]
        Y_base=4*Y_in
        K_total=max(KINV,2*KD)
        add=K_total/M1**2*Icoef*TR*X*X*f
        omega=(Y_base+add)*MPHI*S0/RHOC
        J=JCOEF*TR**3
        rows.append(dict(TR_MeV=TR*1000,X_GeV=X,energy_history_factor=f,Y_base=Y_base,Y_recycle_upper=add,omega_noheat_upper=omega,heat_exponent=J,omega_heat_conservative_upper=omega*math.exp(2*J)))
    ss=np.geomspace(4*M1*M1*(1+1e-9),1e10,1000);maxval=max(direct_sigma_v(s) for s in ss)
    assert maxval<KD
    return dict(KD=KD,KINV=KINV,JCOEF=JCOEF,peak_sampled_sigma_v=maxval,max_relative_energy_derivative_check=max(check),rows=rows)

if __name__=='__main__':print(json.dumps(audit(),indent=2))
