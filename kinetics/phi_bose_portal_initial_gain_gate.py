"""Initial MB inverse-pair portal source and cold scalar gain diagnostic.
No evolving occupation, no actual quantum or cosmological burst.
"""
import math,json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
M=10.;WIDTH=.005305164626087453;A=M*WIDTH;M2=M*M
C=.1**2*1e-5**2*4*math.pi/137.035999084
M1=.1;MCHI2=.02;MPHI=86e-9
GX=2
NODES,WEIGHTS=roots_legendre(24)

def integrated_s(lo,hi,m):
    if hi<=lo:return 0.
    if hi<.5*M2 or lo>1.5*M2:
        s=(hi+lo)/2+(hi-lo)/2*NODES
        return (hi-lo)/2*float(np.dot(WEIGHTS,s*(s+2*m*m)/((s-M2)**2+A*A)))
    x=lo-M2;y=hi-M2
    da=math.atan2(A*(y-x),A*A+x*y)
    val=(y-x)+(M2+m*m)*math.log1p((y*y-x*x)/(x*x+A*A))+(M2*M2+2*m*m*M2-A*A)/A*da
    return val

def source(E,T,m):
    """df/dt per spin and charge at f_dark=0, massless MB e+e-.
    Detailed balance converts inverse into weighted vacuum forward loss.
    """
    if E<=m:return 0.
    p=math.sqrt((E-m)*(E+m));b=(M2-2*m*m)/2
    B=math.sqrt((b-m*m)*(b+m*m))
    resonance_low=(b*b+p*p*m*m)/(E*b+p*B)
    resonance_high=(E*b+p*B)/(m*m)
    pts=[x for x in ((resonance_low-m)/T,(resonance_high-m)/T) if 0<x<80]
    def f(u):
        Ep=m+T*u;pp=math.sqrt((Ep-m)*(Ep+m))
        small=m*m*(E*E+Ep*Ep-m*m)/(E*Ep+p*pp)
        lo=2*m*m+2*small;hi=2*m*m+2*(E*Ep+p*pp)
        return math.exp(-u)*integrated_s(lo,hi,m)
    v=quad(f,0,80,epsabs=0.,epsrel=5e-8,points=pts,limit=180)[0]
    return GX*C*T*math.exp(-(E+m)/T)*v/(192*math.pi**3*E*p)

def scalar_kinematics(k):
    e=math.hypot(k,MPHI)
    est=(M1*M1+MPHI*MPHI-MCHI2*MCHI2)/(2*M1)
    ks=math.sqrt((est-MPHI)*(est+MPHI))
    emin=M1*(est*est+k*k)/(e*est+k*ks)
    emax=M1*(e*est+k*ks)/(MPHI*MPHI)
    return e,emin,emax

def gain_initial(T,alpha):
    k=T*alpha;e,lo,hi=scalar_kinematics(k)
    assert lo>=M1*(1-1e-10) and lo-e>=MCHI2*(1-1e-10)
    cap=min((hi-lo)/T,60.)
    j1=quad(lambda u:source(lo+u*T,T,M1)*T,0,cap,epsabs=0.,epsrel=2e-6,limit=90)[0]
    j2=quad(lambda u:source(lo+u*T-e,T,MCHI2)*T,0,cap,epsabs=0.,epsrel=2e-6,limit=90)[0]
    est=(M1*M1+MPHI*MPHI-MCHI2*MCHI2)/(2*M1)
    ks=math.sqrt((est-MPHI)*(est+MPHI))
    rest=1e-6*((M1+MCHI2)**2-MPHI**2)*ks/(8*math.pi*M1*M1)
    c1=quad(lambda u:source(lo+u*T,T,M1)*(lo+u*T)/(rest*M1)*T,0,cap,epsabs=0.,epsrel=2e-6,limit=90)[0]
    pref=4*rest*M1*M1/(2*ks*k*e)
    gain_bound=pref*.5*c1*c1*max(0.,j1-j2)/(j1*j2)
    delta=max(0.,(j1-j2)/j1)
    t_cut=8*(c1/j1)*delta/(1-4*delta) if delta<.25 else float('inf')
    sharp_bound=min(gain_bound,pref*.5*max(0.,j1-j2)*t_cut*t_cut)
    return {'T_GeV':T,'p_phi_over_T':alpha,'E_parent_min_over_T':lo/T,'int_portal_source_chi1_GeV2':j1,'int_portal_source_chi2_GeV2':j2,'ratio':j1/j2,'normalized_gain':(j1-j2)/(j1+j2),'rest_width_GeV':rest,'mean_parent_residence_GeV_minus1':c1/j1,'integrated_positive_gain_upper_reduced':gain_bound,'markov_positive_duration_upper_GeV_minus1':t_cut,'sharpened_gain_upper_reduced':sharp_bound}

if __name__=='__main__':
    data=[]
    for T in (.1,.2,.4,.6882487043,1.,2.,5.,10.):
        for alpha in (.0001,.0003,.001,.003,.012,.03):
            if T<=.2 and alpha<.003:continue
            data.append(gain_initial(T,alpha))
    print(json.dumps(data,indent=2))
