"""Tree-vacuum 2->4 bound for a single real scalar with V=lam*phi^4/4!."""
import argparse
import json
from pathlib import Path

import numpy as np

from conditional_relic_normalization_reference import conversion
from sm_dark_entropy_exchange_reference import C, MPL
from three_species_decay_onsager_reference import m, mp


def kernel_upper(s,lam,mu=mp):
    """Cauchy and factorized positive four-body phase-space upper bound."""
    x=np.asarray(s)/mu**2
    safe=np.maximum(x,16.)
    t=np.log((safe-1.)/8.)+.125-1./(safe-1.)
    sp=np.log((safe-1.)/3.)**2
    return np.where(x>=16.,lam**4/(mu**2*safe)*(10./24.)*
                    (4*t/(4096*np.pi**5)+6*sp/(2048*np.pi**5)),0.)


def global_upper(lam,mu=mp):
    # log(x/b)^k/x maximized at x=b*exp(k), k=1,2; x>=16.
    t=1/(8*np.e)+1/128
    sp=4/(3*np.e**2)
    return lam**4/mu**2*(10./24.)*(4*t/(4096*np.pi**5)+6*sp/(2048*np.pi**5))


def kinematics_test():
    from scipy.special import roots_legendre
    x,w=roots_legendre(240)
    rng=np.random.default_rng(710)
    checks=[]
    # All energies here are in units of the scalar mass. Using a modest
    # range permits an independent direct subtraction of the two endpoints.
    for _ in range(90):
        A,B=np.exp(rng.uniform(np.log(4.),np.log(20.),2))
        s=(np.sqrt(A)+np.sqrt(B))**2*np.exp(rng.uniform(.0001,1.))
        c=1.+(A+B-s)/2
        h=.5*np.sqrt((1.-4/s)*((s-A-B)**2-4*A*B))
        t1,t2=c+h,c-h
        assert t2<=t1<=0.
        expected=(A-1.)*(B-1.)+(A-B)**2/s
        ar=np.sum(w/(1.-c-h*x)**2)/2
        exact=1./((1.-t1)*(1.-t2))
        checks.append(max(abs(t1*t2/expected-1.),abs(ar/exact-1.)))
    assert max(checks)<2e-11
    x=np.geomspace(16.,1e30,10000)*mp**2
    ratio=np.max(kernel_upper(x,1.)/global_upper(1.))
    assert 0.<ratio<1.
    return dict(max_independent_tu_product_and_angle_relative_defect=max(checks),
                sampled_pointwise_upper_over_closed_form_global_bound=ratio)


def sample_cross_sections():
    """Direct four-body map with all ten tree diagrams (mu=lambda=1)."""
    from itertools import combinations
    from scipy.stats import qmc
    cases=[]
    for S in (18.,24.,40.,100.,500.,2000.):
        estimates=[]
        for seed in (71,907):
            r=qmc.Sobol(8,scramble=True,seed=seed).random_base2(15)
            U=np.sqrt(S)
            a=4.+((U-2.)**2-4.)*r[:,0]
            wB=(U-np.sqrt(a))**2-4.
            b=4.+wB*r[:,1]
            # Three independent isotropic two-body directions.
            def direction(i):
                z=2*r[:,i]-1.
                u=np.sqrt(1-z*z)
                k=2*np.pi*r[:,i+1]
                return np.array([u*np.cos(k),u*np.sin(k),z]).T
            d=direction(2)
            P=np.sqrt(np.maximum((S-a-b)**2-4*a*b,0.))/(2*U)
            Ea=(S+a-b)/(2*U)
            Eb=U-Ea
            def daughters(m2,energy,dp,ix):
                momentum=np.sqrt(m2/4-1.)[:,None]*direction(ix)
                pp=np.sum(momentum*dp,axis=1)
                mass=np.sqrt(m2)
                Er=.5*mass
                def build(sign):
                    en=(energy*Er+sign*P*pp)/mass
                    vec=sign*momentum+(sign*(energy/mass-1.)*pp+P*Er/mass)[:,None]*dp
                    return np.column_stack([en,vec])
                return build(1.),build(-1.)
            A=daughters(a,Ea,d,4)
            B=daughters(b,Eb,-d,6)
            p=np.broadcast_to([U/2,0.,0.,np.sqrt(S/4-1.)],A[0].shape)
            q=np.broadcast_to([U/2,0.,0.,-np.sqrt(S/4-1.)],A[0].shape)
            momenta=[p,q,*(-ev for ev in A+B)]
            inv=lambda k:k[:,0]**2-np.sum(k[:,1:]**2,axis=1)
            phase=((U-2.)**2-4.)*wB*(2*P/U)*np.sqrt((1.-4/a)*(1.-4/b))/(2048*np.pi**5)
            terms=[]
            for i,k in combinations(range(1,6),2):
                Q=momenta[0]+momenta[i]+momenta[k]
                v=inv(Q)
                if i==1:
                    assert v.min()>9.-1e-10
                else:
                    assert v.max()<1e-10
                terms.append(1./(v-1.))
            amp=np.sum(terms,axis=0)
            num=float(np.mean(phase*amp*amp)/(24*S))
            assert 0.<num<kernel_upper(S,1.,1.)
            assert np.max(abs(sum(A+B)-p-q))<2e-11
            assert max(np.max(abs(inv(k)-1.)) for k in A+B)<2e-11
            estimates.append(num)
        cases.append(dict(s_over_mass_squared=S,independent_scramble_estimates_mu_eq_lambda_eq_1=estimates,
                          pointwise_bound_mu_eq_lambda_eq_1=float(kernel_upper(S,1.,1.))))
    return cases


def budget(table,archive,lam,F):
    # Photon floor for a fixed adiabatic visible entropy. We let every
    # represented particle be scalar immediately and remove all 4->2
    # inverse terms, obtaining the maximal growth allowed by this bound.
    t,gs,_=np.loadtxt(table).T
    fee=conversion(table)["Omega_h2_per_GeV_times_yield"]
    T=archive["final_visible_temperature_GeV"]
    phi,chi=archive["conditional_Omega_h2"]
    Y=(phi/mp+chi/m)/fee
    A=np.sqrt(2*C/3)/MPL
    gmax=max(gs[:np.searchsorted(t,T)+1])
    int_bound=(4*C/3)*gmax*T/A
    K=global_upper(lam)
    tau=K*(1.+F)**4*Y*int_bound
    no_stim=K*Y*int_bound
    return dict(archive=archive["archive"],T_GeV=T,quartic_lambda_over_4_factorial=lam,
        final_mode_uniform_F_bound=F,global_vacuum_sigma_v_upper_GeV_minus2=K,
        upper_growth_optical_depth=tau,
        resulting_number_multiplier_upper=1/(1-tau) if tau<1 else None,
        scalar_Omega_h2_upper=fee*mp*Y/(1-tau) if tau<1 else None,
        necessary_uniform_F_cap_for_99_percent_of_reference=max(0.,((1.-fee*mp*Y/(.99*.120))/no_stim)**.25-1.),
        F_at_which_this_comparison_no_longer_gives_a_finite_bound=max(0.,no_stim**(-.25)-1.))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table")
    p.add_argument("normalization")
    p.add_argument("--final-occupation-cap",type=float,default=0.)
    a=p.parse_args()
    assert a.final_occupation_cap>=0.
    d=json.loads(Path(a.normalization).read_text())
    sweep=d["prescribed_vacuum_sweep"]
    hot=next(r for r in sweep["rows"] if r["T_R_GeV"]==.05)
    initial=[dict(archive="all prescribed TR50 source counts at reheating",
             final_visible_temperature_GeV=.05,conditional_Omega_h2=hot["conditional_Omega_h2"]),
             dict(archive="all prescribed TR9.019 source counts at reheating",
             final_visible_temperature_GeV=sweep["vacuum_only_nominal_T_reheat_for_Planck_stable_count_GeV"],
             conditional_Omega_h2=sweep["fine_proxy_at_nominal_T"])]
    results=[budget(a.table,r,L,a.final_occupation_cap)
             for L in (7.65175e-8,1.08515714e-6)
             for r in [*d["at_fixed_archive_counts"],*initial]]
    print(json.dumps(dict(checks=kinematics_test(),sampled_cross_sections=sample_cross_sections(),rows=results,
       scope="Tree diagrams with two quartic vertices only; no sextic contact, Yukawa-loop amplitude, coherent background, or new production. Conserved visible comoving entropy, H at least photon floor. Optional F is a uniform upper bound on all accessible final scalar occupations throughout the interval, not a measured later property."),
       indent=2,allow_nan=False))
