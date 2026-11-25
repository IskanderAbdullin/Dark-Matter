"""Conditional self-blocking / chemistry scales and local burst C22 check.
No claim of actual decay blockade: non-equilibrium medium dispersion,
eta spectral properties, coherence, recoil and gravity must be solved.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad,trapezoid
from late_activation_overlap import H
from direct_activation_consistency import RHO0_CGS,CGS_RHO_TO_EV4
from resonant_mediator_inverse import M,C,HBAR,HBARC,SEC_GYR
from quartic_consistency import lambda_from_sigma
from late_switchon_nr import compute

MPC=3.085677581e24

def track(ps,sigma,rstart=1e8):
    lam=float(lambda_from_sigma(sigma));n0=RHO0_CGS*CGS_RHO_TO_EV4/M;g=lam/(8*M*M)
    delta=2*M*ps/(M+np.sqrt(M*M-2*M*ps))
    ncap=delta/g;r_unlock=(ncap/n0)**(1/3)
    r=np.geomspace(1,rstart,3001);n=n0*r**3;nb=np.minimum(n,ncap);F=nb/n
    Hn=HBAR*H(r-1);mu=g*nb
    loss=2*np.sqrt(3)/(12288*np.pi)*lam**4*nb**4/M**8
    loss_over_H=loss/(n*Hn)
    # Pure coherent weak-coupling tree benchmark; ignore backreaction.
    I=float(trapezoid(loss_over_H,np.log(r)))
    # B-only sound horizon vs fictitious common-velocity locked A+B.
    cs=np.sqrt(mu/M);sel=r>=3400
    sh=float(trapezoid(C/MPC*cs[sel]/H(r[sel]-1),r[sel]))
    sh_lock=float(trapezoid(C/MPC*cs[sel]*np.sqrt(F[sel])/H(r[sel]-1),r[sel]))
    rd=min(r_unlock,rstart)
    dx=C*ps/(M*rd*MPC)*quad(lambda l:np.exp(2*l)/H(np.exp(l)-1),0,np.log(rd),epsrel=1e-8)[0]
    # Maximum redshifted eta radiation relative to one standard nu species
    # at recombination if all parents have decayed by then.
    rho_nu0=7*np.pi**2/120*(4/11)**(4/3)*(2.7255*8.617333262e-5)**4
    dneff_max=(delta/M)*(n0*M)/rho_nu0/1100
    return dict(recoil_eV=ps,sigma_cm2_g=sigma,ncap_eV3=ncap,r_unlock=r_unlock,
                B_fraction_at_r1e8=float(min(1.,(r_unlock/1e8)**3)),
                integrated_coherent_loss_fraction_no_backreaction=I,
                max_loss_per_total_H=float(max(loss_over_H)),max_meanfield_over_m=float(max(mu/M)),
                B_only_contact_sound_horizon_to_equality_Mpc=sh,
                locked_mixture_contact_sound_horizon_to_equality_Mpc=sh_lock,
                kick_displacement_if_conversion_at_unlock_Mpc=dx,
                recombination_delta_neff_upper_if_all_converted=dneff_max,
                scope='Toy instantaneous density cap from NR coherent chemical potential g*nB=delta, not a solved kinematic threshold. Vacuum coherent 4->2 only. Contact-pressure-only sound horizon, not total transfer.')

def local(case,sigma):
    ps=case['recoil_eV'];r=case['r'];p0=case['p0'];lam=float(lambda_from_sigma(sigma));out=[]
    edge=np.asarray(case['edges_B'])
    for key,f in case['B_snapshots'].items():
        F=case['first_crossings'][key]['B_fraction'];t=case['first_crossings'][key]['time_s']
        n=case['N0_eV3']*F;rho=n*M/CGS_RHO_TO_EV4;v=p0*r/M*C
        coeff=rho*sigma*v
        rows=[]
        with np.errstate(divide='ignore',invalid='ignore'):
            for order in [1,2]:
                x=compute(edge,np.asarray(f),1.,order)
                rr=np.array(x['relative_rate_over_n_sigma_vpeak'])*coeff
                rows.append(dict(order=order,relative_rate_per_s=[float(v) if np.isfinite(v) else None for v in rr],
                                 relative_net_over_absolute=x['net_over_absolute']))
        out.append(dict(B_fraction=F,time_s=t,meanfield_over_recoil=lam*n/(8*M*M*ps),
                        meanfield_over_mass=lam*n/(8*M**3),elastic_dilute_rate_per_s=coeff,
                        line_fastest_linear_gain_per_s=case['fastest_linear_growth_s'],
                        orders=rows))
    # Trapezoidal net moment exposure over actual sampled no-C22 trajectory.
    # Missing early intervals are set to zero for this diagnostic.
    tt=np.r_[0.,[x['time_s'] for x in out]]
    rates=np.array([[0.,0.,0.]]+[[(v if v is not None else 0.) for v in x['orders'][-1]['relative_rate_per_s'][:3]] for x in out])
    exposure=trapezoid(abs(rates),tt,axis=0).tolist()
    return dict(recoil_eV=ps,r=r,sigma=sigma,snapshots=out,integrated_abs_log_moment_exposure_screen=exposure,
                scope='Frozen snapshots from no-C22 local burst, not a solution with feedback. Rates are net changes of p4 and hot-number moments, not microscopic collision counts.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bursts',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    result=dict(track=[track(ps,s) for s in [.1,.3,1.25] for ps in [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,.01]],local=[])
    for v in json.loads(Path(a.bursts).read_text()):
        if not v.get('B_snapshots'):continue
        result['local'].append(local(v,.3))
        Path(a.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('LOCAL',v['recoil_eV'],result['local'][-1]['integrated_abs_log_moment_exposure_screen'],flush=True)
    for x in result['track']:
        if x['sigma_cm2_g']==.3:print('TRACK',x['recoil_eV'],x['r_unlock'],x['integrated_coherent_loss_fraction_no_backreaction'],x['B_only_contact_sound_horizon_to_equality_Mpc'],flush=True)
