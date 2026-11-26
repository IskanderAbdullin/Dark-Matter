"""Critical massive Bose physical-DOS IR susceptibility resolution gate.

For hard pair parent P=2GeV the first real-phi energy interval
[mphi,kmin_pair(P,mchi2)] has identically NO on-shell 1<->2 or 2<->2
Yukawa events. It is frozen in that *restricted* collision network.

Its initial physical number can be accurately quadrature-integrated with
18 p-Gauss nodes even as mu_phi approaches mphi, while the physical
entropy-variable Bose susceptibility G_00 can be off by >two orders
of magnitude. An eventual expanding/SM-coupled model cannot substitute
a converged number moment for a converged Jacobian.

BE b(E)=1/expm1((E-mphi+delta_mu)/T), 0<delta_mu, T>0.
G00 = g/(2pi²) int_0^pmax p² dp B0(E)^2 b(1+b),
B0(E)=(1-(E-mphi)/(U-mphi))² is the first Bernstein-P2 hat.
Use stable E-m=p²/(E+m) and asinh(p/sqrt(2m delta)) stretching,
NOT direct sqrt(p²+m²)-m at p much smaller than m.

No condensate delta function, thermal pole, SM bath, FRW or DM relic
is inferred from this reversible excited-state susceptibility.
"""
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from kinetics.phi_bose_physical_cutcell_weak import MPHI,M2
from kinetics.phi_bose_pair_ir_parent_support import global_pair_scalar_support_floor


def first_frozen_cell(P=2.):
    upper=global_pair_scalar_support_floor(P,M2,MPHI)
    assert upper>MPHI
    return float(MPHI),float(upper)


def _scalar_state(p,delta_mu,T,edge):
    """Return physical p² DOS, first Bernstein hat, Bose b, no cancellation."""
    m=MPHI
    E=np.hypot(p,m)
    excitation=p*p/(E+m)       # exact, stable E-m
    eta=(excitation+delta_mu)/T
    b=1/np.expm1(eta)
    B0=(1-excitation/(edge-m))**2
    return p*p/(2*math.pi**2),B0,b


def legendre_physical_moments(delta_mu,n=18,P=2.,T=1.):
    """Actual fixed p-Gauss quadrature used by current Bernstein cell."""
    assert 0<delta_mu<=MPHI and T>0 and n>1
    lo,hi=first_frozen_cell(P)
    pmax=math.sqrt((hi-lo)*(hi+lo))
    z,w=leggauss(n)
    p=(z+1)*pmax/2
    dp=w*pmax/2
    DOS,B,b=_scalar_state(p,delta_mu,T,hi)
    return dict(gauss_order=n,delta_mu_GeV=delta_mu,
                susceptibility_G00_GeV3=float(np.dot(dp,DOS*B*B*b*(1+b))),
                first_number_moment_M0_GeV3=float(np.dot(dp,DOS*B*b)))


def stretched_physical_moments(delta_mu,P=2.,T=1.):
    """Independent adaptive physical quadrature resolving p_core~sqrt(2m delta)."""
    assert 0<delta_mu<=MPHI and T>0
    lo,hi=first_frozen_cell(P)
    pmax=math.sqrt((hi-lo)*(hi+lo))
    core=math.sqrt(2*MPHI*delta_mu)
    umax=math.asinh(pmax/core)
    def integrand(u,index):
        p=core*math.sinh(u)
        dp_du=core*math.cosh(u)
        DOS,B,b=_scalar_state(p,delta_mu,T,hi)
        if index==0:
            return float(dp_du*DOS*B*B*b*(1+b))
        return float(dp_du*DOS*B*b)
    out=[quad(lambda u:integrand(u,i),0,umax,
              epsabs=1e-23,epsrel=2e-11,limit=450)[0]
         for i in (0,1)]
    return dict(delta_mu_GeV=delta_mu,physical_IR_core_momentum_GeV=core,
                transformed_umax=umax,
                susceptibility_G00_GeV3=out[0],
                first_number_moment_M0_GeV3=out[1])



def fixed_geometric_momentum_gauss_grid(P=2.,n_geometric=20,nper=8,
                                        p_min=1e-16):
    """Fixed positive physical p²dp quadrature with log-spaced IR segments.

    DOES NOT depend on Bose chemical potential, so the entropy, moment and
    susceptibility discrete maps can all use the SAME weights; this is
    essential for consistency of implicit quantum KL. The first [0,p_min]
    interval is integrated separately, and all segment GL rules are exact
    for p² capacity (nper>=2).
    """
    lo,hi=first_frozen_cell(P)
    pmax=math.sqrt((hi-lo)*(hi+lo))
    if not (0<p_min<pmax and n_geometric>=4 and nper>=2):
        raise ValueError("Require physical positive geometric quadrature")
    edges=np.r_[0.,np.geomspace(p_min,pmax,n_geometric+1)]
    x,w=leggauss(nper)
    p=(edges[:-1,None]+(edges[1:,None]-edges[:-1,None])*
       (x[None,:]+1)/2).ravel()
    dp=((edges[1:,None]-edges[:-1,None])*w[None,:]/2).ravel()
    assert np.min(p)>0 and np.min(dp)>0
    return p,dp


def composite_physical_moments(delta_mu,P=2.,T=1.,n_geometric=20,nper=8,
                               p_min=1e-16):
    """Positive fixed-grid physical first frozen cell moments and exact DOS capacity."""
    assert 0<delta_mu<=MPHI and T>0
    lo,hi=first_frozen_cell(P)
    pmax=math.sqrt((hi-lo)*(hi+lo))
    p,dp=fixed_geometric_momentum_gauss_grid(P,n_geometric,nper,p_min)
    DOS,B,b=_scalar_state(p,delta_mu,T,hi)
    raw_capacity=float(np.dot(dp,p*p))
    return dict(delta_mu_GeV=delta_mu,nodes=len(p),
                susceptibility_G00_GeV3=float(np.dot(dp,DOS*B*B*b*(1+b))),
                first_number_moment_M0_GeV3=float(np.dot(dp,DOS*B*b)),
                relative_DOS_capacity_error=raw_capacity/(pmax**3/3)-1)

def leading_asymptotic_G00(delta_mu,T=1.):
    """Exact leading δ→0 massive p~sqrt(2mδ) Bose susceptibility singularity."""
    assert 0<delta_mu and T>0
    return T*T*(2*MPHI)**1.5/(8*math.pi*math.sqrt(delta_mu))


def audit():
    lo,hi=first_frozen_cell()
    rows=[]
    for delta in (MPHI,1e-10,1e-12,1e-14,1e-16,1e-18):
        exact=stretched_physical_moments(delta)
        gauss=[legendre_physical_moments(delta,n)
               for n in (18,36,72,144,288,576)]
        rows.append(dict(delta_mu_GeV=delta,adaptive=exact,
                         number_G18_to_reference=gauss[0]["first_number_moment_M0_GeV3"]/
                             exact["first_number_moment_M0_GeV3"],
                         G00_Gauss_to_reference={str(r["gauss_order"]):
                               r["susceptibility_G00_GeV3"]/
                               exact["susceptibility_G00_GeV3"] for r in gauss},
                         G00_IR_asymptote_to_reference=
                               leading_asymptotic_G00(delta)/
                               exact["susceptibility_G00_GeV3"]))
    return dict(schema_version=1,scope="massive real scalar first collision-frozen IR Bernstein-P2 cell physical susceptibility, equilibrium subcritical Bose chemical potential. NOT full five-field finite-time, zero-mode/BEC, SM/FRW or relic",
                first_frozen_phi_energy_edges_GeV=[lo,hi],Ppair_GeV=2.,
                parent_pair_mass_GeV=M2,temperature_GeV=1.,
                rows=rows,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
