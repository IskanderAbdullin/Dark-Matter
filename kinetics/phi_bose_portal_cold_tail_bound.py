"""Low-reheat one-pass cold-scalar count bound, arbitrary vacuum Bose enhancement."""
import math,json
from kinetics.phi_bose_portal_low_reheat_gate import MPL,M,GAMMA,GSTAR,S0,RHOC,MPHI,M1,M2,COUP2
from kinetics.phi_bose_portal_initial_gain_gate import scalar_kinematics
Q_RMS=1.14255e-6 # eV, chosen FD rms diagnostic
T_GAMMA=2.34865e-4 # eV
G_S0=3.91
ALPHA=math.sqrt(2)*Q_RMS/(T_GAMMA*(G_S0/GSTAR)**(1/3))
A=M*GAMMA
sstar=math.hypot(M*M,A)
K=COUP2/(4*math.pi)*(sstar+M*M)/(2*A*A)
COEF=6*K*MPL/(math.pi**4*(2*math.pi**2/45*GSTAR)*1.66*math.sqrt(GSTAR))

def log_y_tail(TR,Emin):
    x=Emin/TR
    return math.log(COEF*TR)+math.log(x+2)-x

def pair_Emin(k,m=M2):
    e=math.hypot(k,MPHI)
    if e>=m:return 2*m
    return 2*m*(m*m+k*k)/(m*e+math.sqrt((m-MPHI)*(m+MPHI))*k)

def logadd(a,b):
    return max(a,b)+math.log1p(math.exp(-abs(a-b)))

def row(TR):
    k=ALPHA*TR;e,Edec,_=scalar_kinematics(k)
    Epair=pair_Emin(k)/2
    assert k<.048 and e<M2 and Edec>=M1
    ld=math.log(2)+log_y_tail(TR,Edec)
    lp=math.log(4)+log_y_tail(TR,Epair)
    log_omega=math.log(MPHI*S0/RHOC)+logadd(ld,lp)
    return dict(TR_MeV=1000*TR,pmax_MeV=1000*k,Edec_min_GeV=Edec,Epair_one_hard_min_GeV=Epair,x_pair=Epair/TR,log10_Y_direct_max=ld/math.log(10),log10_Y_pair_max=lp/math.log(10),log10_omega_cold_upper=log_omega/math.log(10))

if __name__=='__main__':
    from scipy.optimize import brentq
    crossing=brentq(lambda t:row(t)['log10_omega_cold_upper']-math.log10(.06),.025,.05,xtol=1e-12)
    print(json.dumps(dict(alpha=ALPHA,K_e_upper=K,COEF=COEF,rows=[row(t) for t in [.005,.008,.01,.015,.02,.025,.03,.05,.1]],crossing_MeV=1000*crossing),indent=2))
