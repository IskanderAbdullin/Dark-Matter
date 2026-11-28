"""Per-mode quantum decay clocks with held Lorentz-invariant pole masses."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

G=4.
X=.005
KAPPA=.0217
Q=4.1e-6


def intervals(M,m,b,p):
    if M<=m+b:
        return None
    A=(M*M+m*m-b*b)/2
    D=np.sqrt((M*M-(m+b)**2)*(M*M-(m-b)**2))/2
    E=np.hypot(m,p)
    large=A*E+p*D
    # Rationalized lower endpoint; the massless case gives E>=A/(2p)+M²p/(2A).
    lo=(M*M*E*E+D*D)/large
    hi=large/(m*m) if m else np.inf
    return E,lo,hi,2*p*D/(m*m) if m else np.inf,D/M


def parent_integral(l,h):
    w=h-l
    return np.log1p(np.exp(-l)*(-np.expm1(-w))/(1+np.exp(-h)))


def kernel(M,m,b,p,mu=None):
    ranges=intervals(M,m,b,p)
    if ranges is None:
        return None
    E,l,h,dx,P=ranges
    fac=G*M*M/(2*E*p*P)
    if mu is None:
        A=parent_integral(l,l+dx)
        B=A
    else:
        # Wide intervals are truncated in the FD tails, not at a physical pole.
        upper=min(l+dx,max(90.,l+30.))
        A=quad(lambda y:expit(-y)*expit(y-E-mu),l,upper,
               epsabs=1e-14,epsrel=2e-11)[0]
        B=A*(-np.expm1(E+mu))
    yuk=2*P/M*((M+b)**2-m*m)/M**2
    return dict(E=E,lo=l,hi=None if np.isinf(h) else h,
                seed_over_width=fac*A,stimulus_over_width=fac*B,
                vacuum_matched_yukawa_width_ratio=yuk)


def control(M):
    # Integrating over all scalar modes must count one daughter per
    # time-dilated decay of a fixed distribution of parent particles.
    expected=G*M*quad(lambda p:p*p/np.hypot(p,M)*expit(-np.hypot(p,M)),
                       0.,70.,epsabs=1e-12)[0]
    err=[]
    for m,b in ((0.,0.),(5.6e-5,.003),(1e-6,.0048),(.0019,.003)):
        f=lambda l:np.exp(3*l)*kernel(M,m,b,np.exp(l))['seed_over_width']
        q=quad(f,np.log(1e-12),np.log(80.),epsabs=1e-11,epsrel=5e-9)[0]
        err.append(abs(q/expected-1))
    assert max(err)<1e-7
    dr=[]
    for m,b,p in ((0.,0.,4.1e-6),(.00005,.003,4.1e-6),(.0019,.003,.01)):
        E,l,h,dx,P=intervals(M,m,b,p)
        mu=-.55*E
        upper=min(h,90.)
        dst=quad(lambda y:expit(-y)-expit(-y+E+mu),l,upper,
                 epsabs=1e-13,epsrel=5e-10)[0]
        q=kernel(M,m,b,p,mu)
        dr.append(abs(q['stimulus_over_width']/(G*M*M*dst/(2*E*p*P))-1))
    assert max(dr)<3e-7
    df=kernel(M,5.6e-5,0.,1e-10)['seed_over_width']
    exp=G*(M/5.6e-5)**3*expit(-(M*M+5.6e-5**2)/(2*5.6e-5))
    assert abs(df/exp-1)<5e-10
    return dict(total_seed_time_dilated_count_max_defect=max(err),
                independent_FD_gain_loss_max_defect=max(dr),
                finite_mass_zero_momentum_max_defect=abs(df/exp-1))


def run():
    mlist=(0.,1e-9,2e-6,4.1e-6,1e-5,np.sqrt(7.65175e-8/24),
           .00015,.001,.002,.0049)
    values=[]
    annual=KAPPA*X*X
    for p in (1e-6,Q,1e-5):
        ref=kernel(X,0.,0.,p)['seed_over_width']
        for b in (0.,.003,.0048,.0051):
            for m in mlist:
                e=kernel(X,m,b,p)
                if e is None:
                    values.append(dict(momentum_over_T=p,scalar_mass_over_T=m,fermion_mass_over_T=b,
                                       channel_open=False))
                    continue
                occ=kernel(X,m,b,p,mu=-e['E']/2)
                under=kernel(X,m,b,p,mu=-2*e['E'])
                values.append(dict(momentum_over_T=p,scalar_mass_over_T=m,fermion_mass_over_T=b,
                                   channel_open=True,parent_energy_lo_over_T=e['lo'],parent_energy_hi_over_T=e['hi'],
                                   empty_ST_rate_per_H_fixed_width=e['seed_over_width']*annual,
                                   empty_rate_relative_to_massless_daughters_fixed_width=e['seed_over_width']/ref,
                                   illustrative_yukawa_width_ratio=e['vacuum_matched_yukawa_width_ratio'],
                                   empty_rate_per_H_yukawa_width=e['seed_over_width']*annual*e['vacuum_matched_yukawa_width_ratio'],
                                   half_scalar_energy_fermion_negative_mu_stimulus_per_H=occ['stimulus_over_width']*annual,
                                   twice_scalar_energy_fermion_negative_mu_stimulus_per_H=under['stimulus_over_width']*annual))
    return dict(scope="Instantaneous narrow on-shell benchmark, pole energies E²=p²+m² with unit residues, "
                      "fixed FD parent and empty or held FD daughter. No in-medium widths, different "
                      "dispersion slopes, chemical evolution, FRW solution or actual scalar thermal mass.",
                M_over_T=X,vacuum_width_over_H_pivot=annual,daughter_mode_pivot=Q,validation=control(X),cases=values)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
