"""Identify the source of the frozen Hartree-mass difference, without an evolution."""
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from finite_mass_qed_cold_sensitivity_pilot import GRIDS,K,load_all
from finite_mass_qed_timelike_pole_angular_reference import angle


LAM=7.65175e-8


def quadrature(edges,f,order):
    x,y=roots_legendre(order)
    l,r=edges[:-1,None],edges[1:,None]
    q=(l+r)/2+(r-l)*x/2
    return q,(r-l)*y/2*q*q*f[:,None]/(2*np.pi**2)


def reconstruct(arc,native,match_count):
    e=np.array(arc['scalar_momentum_edges'])
    midpoint=np.array(arc['scalar_momentum_nodes'])
    f=np.array(arc['endpoint_scalar_occupations'])
    mass0=arc['x_end']*1e-9
    accepted=np.array(arc['endpoint_scalar_signed_shell_number_over_T3'])
    cells=np.diff(e**3)/(6*np.pi**2)
    if not np.array_equal(f*cells,accepted):
        raise ArithmeticError('Raw shell occupation and counts differ')
    dense=lambda m,n=accepted: np.sum(n/np.hypot(m,midpoint))
    n0=sum(accepted)
    if abs(native*native-mass0*mass0-.5*LAM*dense(native))/native**2>3e-13:
        raise ArithmeticError('Saved Hartree gap equation not recovered')
    masses=[]
    for order in (24,40):
        q,W=quadrature(e,f,order)
        if abs(W.sum()-n0)>2e-14*n0:
            raise ArithmeticError('Exact shell quadrature changed number')
        def solve(scale):
            return brentq(lambda m:m*m-mass0*mass0-.5*LAM*np.sum(scale*W/np.hypot(m,q)),
                          native*.85,native*1.15,xtol=1e-20)
        masses.append([solve(1.),solve(match_count/n0)])
    prev=(np.array(masses[1])-masses[0])/masses[1]
    if np.max(abs(prev))>2e-13:
        raise ArithmeticError('Smooth shell mass integral not converged')
    m,common=masses[-1]
    E=np.hypot(m,q)
    I=np.sum(W/E)
    deriv=(.5*LAM*I)/(m*(2*m+.5*LAM*np.sum(W*m/E**3)))
    h=1e-4
    check=(np.log(solve(1+h))-np.log(solve(1-h)))/(np.log1p(h)-np.log1p(-h))
    if abs(deriv-check)>1e-8:
        raise ArithmeticError('Implicit density derivative did not close')
    # Evaluate the tail fraction without straddling a histogram edge.
    all_e=np.unique(np.r_[e,1e-5])
    q1,W1=quadrature(all_e,f[np.searchsorted(e,(all_e[:-1]+all_e[1:])/2)-1],40)
    frac=np.sum((W1/np.hypot(m,q1))[q1<1e-5])/I
    return dict(archived_mass=native,scalar_total_number_over_T3=n0,extra_vacuum_mass=mass0,
                exact_step_gap_mass=m,common_count_step_gap_mass=common,
                exact_step_minus_native_relative=(m-native)/native,
                working_density_elasticity=deriv,tadpole_fraction_q_below_1e_5=frac)


def run():
    F=load_all()
    want=sum(F['third'][0]['endpoint_scalar_signed_shell_number_over_T3'])
    cases={g:reconstruct(F[g][0],F[g][2],want) for g in GRIDS}
    a,b=[cases[g] for g in GRIDS]
    nr=a['scalar_total_number_over_T3']/b['scalar_total_number_over_T3']
    dm=a['exact_step_gap_mass']/b['exact_step_gap_mass']-1
    resid=a['common_count_step_gap_mass']/b['common_count_step_gap_mass']-1
    if abs(dm-(nr**(1/3)-1))>2e-6 or abs(resid)>2e-6:
        raise ArithmeticError('Equal count did not control the archived mass drift')
    _,ch,m,w,R=F['fine']
    rows=[]
    for q in (3.5e-6,8.2e-6):
        def coeff(m0):
            return angle(ch,K,w,R,np.hypot(m0,q),'T',6,(6,12),q)[0]
        C=coeff(m)
        step=1e-6*m
        z=[]
        for h in (5*step,step):
            z.append((np.log(coeff(m+h))-np.log(coeff(m-h)))/(np.log(m+h)-np.log(m-h)))
        if np.max(abs(z[0]/z[1]-1))>4e-6:
            raise ArithmeticError('Local dispersion response unresolved')
        slope=z[-1]
        neutral_tol=np.log1p(.01)/abs(slope[0])*1./a['working_density_elasticity']
        rows.append(dict(q=q,raw_A_B_D=C.tolist(),d_log_coeff_d_log_mass=slope.tolist(),
                         local_fractional_density_tolerance_for_1_percent_raw_forward=neutral_tol))
    times=[]
    for a0,b0 in zip(F['fine'][0]['checkpoints'],F['third'][0]['checkpoints'],strict=True):
        if a0['x']!=b0['x'] or min(a0['scalar_yield'],b0['scalar_yield'])<=0:
            raise ArithmeticError('Archived checkpoints are not paired')
        times.append(dict(x=a0['x'],fine_to_third_yield=a0['scalar_yield']/b0['scalar_yield'],
                          fine_to_third_pole=a0['scalar_pole_over_T']/b0['scalar_pole_over_T']))
    if abs(times[-1]['fine_to_third_yield']-nr)>4e-14:
        raise ArithmeticError('Endpoint number and archived yield ratios disagree')
    return dict(scope='Frozen Hartree bookkeeping and local conditional sensitivity, not an evolved remesh.',
                common_count_over_T3=want,counts_ratio=nr,relative_step_mass_difference=dm,
                common_count_relative_step_mass_difference=resid,cases=cases,archive_checkpoint_ratios=times,local_response=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
