"""Number-vs-momentum envelope for conditional collisionless 86 eV scalar states."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import mp
from warm_prompt_reheating_expanding_slow_reference import visible


def packed(bg,ai,n0,cap,rad,order):
    # n0 is comoving with a_today=1, for one real scalar degree of freedom.
    q=(6*np.pi*np.pi*n0/cap)**(1/3)
    x,w=roots_legendre(rad)
    z=(1+x)/2
    ww=3*z*z*w/2
    assert abs(ww.sum()-1.)<3e-14
    result=bg.streams(q*z,mp,ai,bg.aeq,order)
    assert np.all(np.diff(result)>0)
    return dict(cap=cap,q_edge_over_Tgamma0=q/bg.T0,
                q_median_over_Tgamma0=2**(-1/3)*q/bg.T0,
                rms_q_over_Tgamma0=np.sqrt(3/5)*q/bg.T0,
                maximal_count_fraction_below_q_over_Tgamma0_point01=min(1.,(.01*bg.T0/q)**3),
                maximal_count_fraction_below_q_over_Tgamma0_point1=min(1.,(.1*bg.T0/q)**3),
                minimum_mean_path_Mpc=float(ww@result))


def run(eos,path):
    with gzip.open(path,'rt') as f:
        path=json.load(f)['arms'][-1]['trajectory']
    st=path[-1]
    T=st['T_visible']
    bg=Background(eos)
    ai=bg.scale(T)
    _,s=visible(bg,T)
    conv=conversion(eos)['Omega_h2_per_GeV_times_yield']
    entcom=s*ai**3
    direct=(2*np.pi*np.pi/45)*bg.gs0*bg.T0**3
    assert abs(entcom/direct-1)<2e-14
    rows=[]
    for omega in (.06,.12):
        Y=omega/(conv*mp)
        n=entcom*Y
        opts=[]
        for F in (.52,1.,10.,100.,1e3,1e6,1e9):
            a=packed(bg,ai,n,F,18,5)
            b=packed(bg,ai,n,F,40,8)
            assert abs(a['minimum_mean_path_Mpc']/b['minimum_mean_path_Mpc']-1)<3e-6
            opts.append(b)
        specials=[]
        for path_Mpc in (1.,.3,.1,.03):
            z=brentq(lambda z:packed(bg,ai,n,np.exp(z),40,8)['minimum_mean_path_Mpc']-path_Mpc,
                     -15.,30.,xtol=1e-11)
            specials.append(dict(proposed_mean_path_upper_Mpc=path_Mpc,
                                 necessary_occupation_cap=np.exp(z)))
        fcap1=packed(bg,ai,n,1.,40,8)
        for qstar in (.01,.1):
            specials.append(dict(proposed_median_q_over_Tgamma0_upper=qstar,
                                 necessary_occupation_cap=.5*(fcap1['q_edge_over_Tgamma0']/qstar)**3))
        rows.append(dict(target_omega=omega,target_Y=Y,cap_scan=opts,illustrative_path_or_median_targets=specials))
    # Thermal shape control for reference only, not the minimum spectrum.
    Td=st['light_scalar_T'][1]
    peak=st['phase_factors'][1]
    g,w=roots_legendre(80)
    z=(g+1)*35.
    p=Td*z
    E=np.hypot(mp,p)
    ph=peak*np.exp(-(p*p/(E+mp))/Td)
    W=35*Td*w*p*p*ph/(1-ph)/(2*np.pi*np.pi)
    assert abs(W.sum()/(st['yields'][1]*s)-1)<2e-7
    path_mean=(W@bg.streams(p*ai,mp,ai,bg.aeq,8))/W.sum()
    return dict(scope="Fixed-background collisionless geometric phase-space packing for one real scalar, "
                      "not a transfer prediction, cutoff requirement, cosmological likelihood or mechanism. "
                      "Path and median targets are illustrative, not observational limits. "
                      "No bound on actual cap, coherent component, or collisional cooling.",
                visible_T_GeV=T, reference_equality_a=bg.aeq,
                fixed_targets=rows, modeled_thermal_shape_control=dict(
                    fitted_rest_occupation=peak/(1-peak),
                    conditional_omega_if_counts_freeze=st['yields'][1]*mp*conv,
                    fixed_background_mean_path_Mpc=float(path_mean)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('trajectory',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.trajectory),indent=2,allow_nan=False))
