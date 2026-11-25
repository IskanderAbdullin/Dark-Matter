"""Independent kinematic and weak-coupling checks of derivative production.

Conventions: L=(d_mu A/F) bar(N2) gamma^mu (cV+cA gamma5) N1+h.c.
One Dirac transition, spin averaged. F=2 f/C in existing pseudoscalar
rate_mapping. Shift symmetry is a property of this EFT, not of a lone
on-shell-equivalent Yukawa coupling. No claim of UV completion.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import gamma,jv
from scipy.integrate import solve_ivp
from yukawa_rate_bath import rate_mapping

MP=1.220890e19
MPL=2.435e18
M=86000.
m=86e-9
GSTAR=106.75
R=0.9
K=.00685
RHO_OVER_S_DM=.12*1.053672e-5/2891.2 # GeV

def widths(M,r,m,F):
    m2=r*M
    p=np.sqrt(((M-m2)**2-m*m)*((M+m2)**2-m*m))/(2*M)
    AV=(M-m2)**2*((M+m2)**2-m*m)/F**2
    AA=(M+m2)**2*((M-m2)**2-m*m)/F**2
    return p*AV/(8*np.pi*M*M),p*AA/(8*np.pi*M*M)

def traces(M,r,m,F):
    """Independent numerical gamma-matrix traces in the parent rest frame."""
    I=np.eye(2);Z=np.zeros((2,2))
    pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    g=[np.block([[I,Z],[Z,-I]])]+[np.block([[Z,s],[-s,Z]]) for s in pauli]
    g5=1j*g[0]@g[1]@g[2]@g[3];I4=np.eye(4)
    slash=lambda p:sum(v*x for v,x in zip(np.asarray(p)*[1,-1,-1,-1],g))
    m2=M*r;p=np.sqrt(((M-m2)**2-m*m)*((M+m2)**2-m*m))/(2*M)
    E2=np.hypot(m2,p);EA=np.hypot(m,p)
    p1=[M,0,0,0];p2=[E2,0,0,-p];k=[EA,0,0,p]
    out=[]
    for current in [I4,g5]:
        a=slash(k)@current/F;abar=g[0]@a.conj().T@g[0]
        mat=(slash(p2)+m2*I4)@a@(slash(p1)+M*I4)@abar
        amp=np.trace(mat).real/2
        out.append(p*amp/(8*np.pi*M*M))
    return out

def harmonic_yield(F,theta,m=m,g=GSTAR):
    # x=mt, radiation H=1/(2t), T^2=m MP/(3.32 sqrt(g) x).
    # theta=theta_i*2^(1/4)*Gamma(5/4)*x^(-1/4)*J_(1/4)(x).
    C=2**.25*gamma(1.25)*np.sqrt(2/np.pi)
    Tscale=np.sqrt(m*MP/(3.32*np.sqrt(g)))
    sscale=2*np.pi**2*g*Tscale**3/45
    return .5*m*m*F*F*theta*theta*C*C/sscale

def cosine_yield(F,theta,xend=1000.,rtol=2e-9):
    x0=1e-4
    init=[theta-np.sin(theta)*x0*x0/5,-2*np.sin(theta)*x0/5]
    fun=lambda x,y:[y[1],-1.5*y[1]/x-np.sin(y[0])]
    t=np.linspace(xend-50.,xend,1001)
    sol=solve_ivp(fun,(x0,xend),init,t_eval=t,rtol=rtol,atol=rtol*1e-3,max_step=.5)
    assert sol.success
    a,v=sol.y
    # Stable 1-cos, late period average of comoving energy.
    q=t**1.5*(.5*v*v+2*np.sin(a/2)**2)
    Tscale=np.sqrt(m*MP/(3.32*np.sqrt(GSTAR)));sscale=2*np.pi**2*GSTAR*Tscale**3/45
    y=float(q.mean()*m*m*F*F/sscale)
    return dict(theta_i=theta,energy_per_entropy_GeV=y,fraction_DM=y/RHO_OVER_S_DM,
        enhancement_over_harmonic=y/harmonic_yield(F,theta),
        late_comoving_energy_relative_span=float((q.max()-q.min())/q.mean()),xend=xend,rtol=rtol)

def screen():
    p=rate_mapping(M,R,K,.03,47.,GSTAR)
    F=2*p['f_over_C_GeV']
    H=p['H_ref_GeV'];wa=widths(M,R,m,F)
    # Test finite daughter-mass Ward relation well away from threshold;
    # the actual tiny mA difference is beneath float precision.
    trace_checks=[]
    for ma in [0.,.01*M,.05*M,.099*M]:
        analytic=widths(M,R,ma,F);num=traces(M,R,ma,F)
        trace_checks.append(dict(mA_over_M=ma/M,analytic=list(analytic),trace=num,
           maximum_relative_error=float(np.max(abs(np.array(num)/analytic-1)))))
    ys=[cosine_yield(F,x) for x in [.01,.1,.5,1.,2.,3.,np.pi-.001]]
    # Generic photon anomaly with c_gamma=1: illustrative, not a prediction.
    alpha=1/137.035999084
    gag=alpha/(2*np.pi*F);gagwidth=gag*gag*m**3/(64*np.pi)
    return dict(mapping=p,F_GeV=F,widths_GeV=list(wa),width_over_H_ref=wa[1]/H,
       trace_checks=trace_checks,pNGB_quartic=-m*m/F**2,
       harmonic_mass_over_H_equals_1_T_GeV=float(np.sqrt(m*MP/(1.66*np.sqrt(GSTAR)))),
       T_at_mt1_GeV=float(np.sqrt(m*MP/(3.32*np.sqrt(GSTAR)))),
       required_DM_energy_per_entropy_GeV=RHO_OVER_S_DM,misalignment=ys,
       theta_small_angle_for_1percent_DM=float(np.sqrt(.01*RHO_OVER_S_DM/harmonic_yield(F,1))),
       theta_small_angle_for_100percent_DM=float(np.sqrt(RHO_OVER_S_DM/harmonic_yield(F,1))),
       illustrative_photon_lifetime_s=6.582119569e-25/gagwidth,
       scope='Exact vacuum Ward kinematics and homogeneous radiation-era cosine misalignment only; finite-temperature kinetic dispersion, UV completions, symmetry-breaking mixing and initial conditions remain open.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=screen()
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps(r,indent=2))
