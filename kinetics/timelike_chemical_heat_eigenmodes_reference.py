"""Local coupled chemical/heat eigenmodes for a restricted MB portal."""
import argparse,json
import numpy as np
from scipy.special import roots_legendre,kve
from timelike_portal_onsager_coverage_reference import integrated,mchi,me,moment_kernel
from sm_dark_entropy_exchange_reference import StandardModel,C,MPL

leptons=(me,.1056583755,1.77686)
P=np.diag([2.,1.])


def rates(T):
    return integrated(T,2*mchi,max(60*T,10+30*T),
                      kernel=lambda E,t:sum(moment_kernel(E,t,m) for m in leptons))


def thermo(T):
    z,w=roots_legendre(120)
    v=np.sqrt(60.)*(z+1)/2
    q=v*v
    E=T*np.hypot(q,mchi/T)
    wt=np.sqrt(60.)*v*w*q*q*4/(2*np.pi**2)*T**3*np.exp(-E/T)
    A,B,D=(wt.sum(),wt@(E/T),wt@((E/T)**2))
    analytic=4*mchi*mchi*T*kve(2,mchi/T)*np.exp(-mchi/T)/(2*np.pi**2)
    assert abs(A/analytic-1)<1e-10
    # Dark scalar is massless and thermally equilibrates but
    # does not contribute to chemical-number susceptibility.
    return np.array([[A,B],[B,D+4*C*T**3]]),wt@E+C*T**4


def trial(T,eos):
    rho,_,_,cs=eos.state(T)
    chi,rd=thermo(T)
    H=np.sqrt((rho+rd)/3)/MPL
    M=rates(T)
    L=np.array([[M[0],M[1]/T],[M[1]/T,M[2]/(T*T)]])
    Sd=np.linalg.inv(chi)
    S=Sd+np.diag([0.,1/cs])
    jac=P@L@P@S
    # Similar to a symmetric PSD matrix via the entropy metric.
    sq=np.linalg.cholesky(S)
    sym=sq.T@P@L@P@sq
    ev=np.linalg.eigvalsh(sym)
    assert np.all(ev>0)
    assert np.max(abs(np.sort(np.linalg.eigvals(jac).real)/ev-1))<3e-9
    # Direct finite difference of the nonlinear MB rate, including both
    # energy and charge forcing, around the common-T reference.
    base=chi[0,0]
    def rhs(y):
        eta,dt=Sd@y
        sm=-y[1]/cs
        b=rates(T*(1+dt))*np.exp(2*eta)-rates(T*(1+sm))
        return -np.array([2*b[0],b[1]/T])
    h=base*1e-5
    jacref=np.stack([(rhs(-v)-rhs(v))/(2*h) for v in (h*np.eye(2))],axis=1)
    check=float(np.linalg.norm(jacref-jac)/np.linalg.norm(jac))
    assert check<5e-7
    return dict(T_GeV=T,eigenmodes_per_H=(ev/H).tolist(),
                pair_forward_per_particle_over_H=float(2*M[0]/(base*H)),
                physical_dimensionless_net_N_Q_correlation=float(L[0,1]/np.sqrt(L[0,0]*L[1,1])),
                rate_matrix_condition=float(ev[-1]/ev[0]),
                crosscheck_finite_difference_fraction=check)


def run(table):
    eos=StandardModel(table)
    rows=[trial(t,eos) for t in (.07,.15,.3,.45,.55,.6882487043,.8,1.)]
    return dict(rows=rows,scope=("Local thermal linear response, one symmetric MB stable species plus thermal "
                    "relativistic scalar, SM EOS, timelike e/mu/tau only. No spacelike heat, hadrons, "
                    "Pauli 2-to-2, other dark or real mediator dynamics, or relic result."))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table")
    print(json.dumps(run(p.parse_args().table),indent=2,allow_nan=False))
