"""MB timelike number/heat moment response with explicit spectral coverage.

References: CLEO III PRD 76 072008 (2007) and upstream draft PR #1
phi_bose_vector_portal_annihilation_resonance.py. No continuous hadronic
R(s) is assumed or inferred. The unity spectral kernel below is a coverage
weight only, not a hadronic prediction.
"""
import json
import math
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import kve


V=10.
width=.005305164626087453
me=.00051099895
mchi=.02
mpi=.13957039
alpha=1/137.035999084
gprod=.1*1e-5*np.sqrt(4*np.pi*alpha)
factor=gprod**2/(12*np.pi)*4/(32*np.pi**4)


def moment_kernel(E,T,lf=me):
    s=E*E
    sm=4*mchi*mchi
    sf=4*lf*lf
    if not s>max(sm,sf):
        return np.zeros(3)
    # Numerator of sigma times invariant flux; BW denominator is separate.
    a=factor*(s+2*mchi*mchi)*E*E*E*np.sqrt((s-sm)/s)
    if lf:
        a*=(1+2*lf*lf/s)*np.sqrt((s-sf)/s)
    x=E/T
    K0,K1=kve(0,x)*np.exp(-x), kve(1,x)*np.exp(-x)
    return a*np.array([T*K1,2*T*T*K1+T*E*K0,
                      (6*T**3+T*E*E)*K1+3*T*T*E*K0])


def coverage_kernel(E,T):
    # Heaviside on 2 pion threshold ONLY. Does not assume actual R=1;
    # it weights a trial constant multiplier to assess data coverage.
    return moment_kernel(E,T,0.) if E>2*mpi else np.zeros(3)


def integrated(T,a,b,w=width,kernel=moment_kernel,coord="log"):
    if b<=a: return np.zeros(3)
    bw=V*w
    p=[a,b,V,2*mchi,2*mpi]
    p += [V+j*w for j in (-300,-100,-30,-10,-3,-1,1,3,10,30,100,300)]
    p += [T*j for j in (.1,.2,.5,1,2,3,5,7,10,15,22,30,50)]
    x=sorted(set(y for y in p if a<=y<=b))
    if coord=="log":
        def fn(t):
            E=np.exp(t)
            return kernel(E,T)*2*E*E/((E*E-V*V)**2+bw*bw)
        x=np.log(x)
    else:
        def fn(th):
            E=np.sqrt(max(V*V+bw*np.tan(th),0.))
            return kernel(E,T)/bw
        x=np.arctan((np.array(x)**2-V*V)/bw)
    return sum((quad_vec(fn,l,h,epsrel=2e-10,epsabs=1e-180,limit=200)[0]
                for l,h in zip(x[:-1],x[1:])),np.zeros(3))


def leptonic_row(T):
    lo=2*mchi
    hi=max(60*T,V+30*T)
    p=integrated(T,lo,hi)
    ref=integrated(T,lo,hi,coord="theta")
    err=float(np.max(abs(p/ref-1)))
    assert err<2e-7
    narrow=np.pi/(V*width)*moment_kernel(V,T)
    cs=V*width
    exact_pole=integrated(T,V-.01,V+.01,coord="theta")
    windows={"near_pole_fraction":(exact_pole/p).tolist(),
             "signed_delta_subtraction_remainder_fraction":((p-narrow)/p).tolist(),
             "finite_window_plus_shoulder_closure_fraction":float(np.max(abs(
                    (exact_pole+integrated(T,lo,V-.01)+integrated(T,V+.01,hi))/p-1)))}
    var=p[2]/p[0]-(p[1]/p[0])**2
    normalized_determinant=1-p[1]**2/(p[0]*p[2])
    assert var>0 and windows["finite_window_plus_shoulder_closure_fraction"]<2e-8
    beta=1/T
    eps=beta*2e-5
    prev=integrated(1/(beta-eps),lo,hi)
    nxt=integrated(1/(beta+eps),lo,hi)
    d_beta=-(nxt[:2]-prev[:2])/(2*eps)
    derivative=float(np.max(abs(d_beta/p[1:]-1)))
    assert derivative<5e-6
    # Force vector (delta alpha = sum(mu/T), delta inverse temperature).
    frac=.0001
    d=np.array([.6*frac,.4*frac/T]) # evaluated dark minus SM, coord -beta
    f=integrated(1/(beta-d[1]/2),lo,hi)*np.exp(d[0]/2)
    r=integrated(1/(beta+d[1]/2),lo,hi)*np.exp(-d[0]/2)
    exact=f[:2]-r[:2]
    pred=np.array([[p[0],p[1]],[p[1],p[2]]])@d
    linear=float(np.max(abs(exact/pred-1)))
    assert linear<2e-6 and d@exact>=0
    return dict(T=T,N_Q_Q2=p.tolist(),nwa_over_bw=(narrow/p).tolist(),**windows,
                coordinate_relative_error=err,beta_derivative_relative_error=derivative,
                linear_response_relative_error=linear,energy_variance_GeV2=var,
                normalized_onsager_determinant=float(normalized_determinant))


def scoped_hadron_coverage(T):
    lo=2*mpi
    hi=max(60*T,V+30*T)
    args=dict(kernel=coverage_kernel)
    all_=integrated(T,lo,hi,**args)
    cov=integrated(T,6.964,10.538,**args)
    point=integrated(T,9.986,10.006,**args)
    under=integrated(T,lo,6.964,**args)
    over=integrated(T,10.538,hi,**args)
    assert np.max(abs((under+cov+over)/all_-1))<3e-9
    max_positive=integrated(T,9.99,10.01,**args)
    return dict(T=T,below_first_CLEO_energy=(under/all_).tolist(),
                across_CLEO_energy_envelope=(cov/all_).tolist(),
                local_point_proxy_window=(point/all_).tolist(),
                symmetric_pole_window=(max_positive/all_).tolist(),
                above_last_CLEO_energy=(over/all_).tolist())


def main():
    T=(.07,.15,.3,.45,.55,.6882487043,.8,1.,1.5)
    return dict(leptonic=[leptonic_row(x) for x in T],
                hadron_coverage=[scoped_hadron_coverage(x) for x in T],
                roles={"leptonic":"vacuum MB rates with physical electron mass, not the SM bath",
                       "hadronic":"positive response kernel coverage, not measured hadron rates"},
                cleo_DOI="https://doi.org/10.1103/PhysRevD.76.072008",
                RIS_primary="https://arxiv.org/abs/2309.16615")
if __name__=="__main__":
    print(json.dumps(main(),indent=2,allow_nan=False))
