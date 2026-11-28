"""Global energy-weighted ceilings on the specified vacuum recycling kernels."""
import json
import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude
from archived_heavy_recycling_gate_reference import angular_integral
from three_species_decay_onsager_reference import m, M, mp, y
from warm_same_charge_heavy_recycling_reference import exact_angular


def geom_sq(s):
    return m**4 + m*m*(s-2*M*M-2*mp*mp) + (M*M-mp*mp)**2


def finite_scalar_cells(n_s, n_x, end):
    """Bound every s slab by an outward-padded positive angle Riemann sum."""
    s = np.geomspace(4*M*M, end, n_s+1)
    # Exact x coverage is [0,1]; reflect since the amplitude is even.
    x = np.linspace(0.,1.,n_x+1)
    zlo,zhi = (x[:-1]**2)[None,:],(x[1:]**2)[None,:]
    zarg=np.clip(.5,zlo,zhi)
    perp=zarg*(1-zarg)
    l,h=s[:-1,None],s[1:,None]
    alo,ahi=m*m+l/2-M*M-mp*mp,m*m+h/2-M*M-mp*mp
    plo,phi=np.maximum(l/4-M*M,0.),h/4-M*M
    klo,khi=l/4-mp*mp,h/4-mp*mp
    # L=(M+m)*A-2*M*k^2*z grows in s and falls in z.
    Lmin=(M+m)*alo-2*M*klo*zhi
    Lmax=(M+m)*ahi-2*M*khi*zlo
    boundL=np.maximum(Lmin**2,Lmax**2)
    D2_min=4*plo*klo
    denom=geom_sq(l)+D2_min*(1-zhi)
    boundA=8*phi*y**4*(h*khi**2*perp+boundL)/(denom**2)
    # Four-factor padding safely dominates observed IEEE evaluation
    # errors for these positive finite slab sums; no derivative guesses.
    partial=(boundA @ (2*np.diff(x)))*(1+1e-11)
    coeff=M*M/np.pi*np.sqrt(1-4*M*M/h[:,0])/l[:,0]**2 * partial
    return float(max(coeff)), [float(s[np.argmax(coeff)]),float(s[np.argmax(coeff)+1])]


def scalar_tail(end):
    # Apply the analytic geometric virtual denominator and beta<=1.
    return 4*M*M/end * y**4/(np.pi*np.sqrt(geom_sq(end))) * (1+1e-11)


def stable_all():
    d=M*M-m*m
    B=(M+m)**2-mp*mp
    I=2*(1+B/d)**2
    # Apply K_o(s) and K_s(s) to energy products, with the identical
    # incoming factor kept separately. Low invariant masses have
    # nonnegative crossed spin product, hence K_s <= K_o there.
    a,b=4*M*M,.095
    c0=lambda s:-s*s/8+s*(m*m+M*M+m*M)-2*m*m*M*M
    assert min(c0(a),c0(b))>0. and b>=5*M*M
    C=y**4*M*M/(4*np.pi)*I*(1+1e-11)
    low=C/(25*M**4*np.sqrt(5))
    # Beyond b use the independent Cauchy bound K_s <= 2 K_o. The
    # beta/s^2 envelope is decreasing above 5 M^2.
    high=2*C*np.sqrt(1-4*M*M/b)/b**2
    return max(low,high), low, high


def run():
    controls=[]
    for q in (170,340,680):
        z,where=finite_scalar_cells(q,3*q//2,.8)
        controls.append(dict(s_cells=q,angle_cells=3*q//2,
                             scalar_finite_upper=z,peak_slab=where))
    tail=scalar_tail(.8)
    test,k_low,k_high=stable_all()
    C=max(tail, controls[-1]["scalar_finite_upper"], test)
    s=np.geomspace(4*M*M*(1+1e-8),16.,1800)
    x,w=roots_legendre(220)
    direct=(amplitude(s[:,None],x,light=M,heavy=m,scalar=mp)@w)*(
        np.sqrt(1-4*M*M/s)*M*M/(np.pi*s*s))
    factor=y**4*M*M/(4*np.pi*s*s)*np.sqrt(1-4*M*M/s)
    stable=angular_integral(s)*factor
    same=exact_angular(s)*factor
    assert max(direct)<max(tail,controls[-1]["scalar_finite_upper"])
    assert max(stable)<test and max(same)<test
    assert test<controls[-1]["scalar_finite_upper"] and tail<C
    return dict(simultaneous_energy_density_coefficient_GeVminus2=C,
                stable_all_charge_coefficient_GeVminus2=test,
                stable_low_and_high_upper_GeVminus2=[k_low,k_high],
                scalar_analytic_tail_GeVminus2=tail,
                scalar_finite_cells=controls,
                sampled_scalar_peak_GeVminus2=float(max(direct)),
                sampled_stable_opposite_peak_GeVminus2=float(max(stable)),
                sampled_stable_same_peak_GeVminus2=float(max(same)),
                scope="Global vacuum rate ceilings for scalar or stable-pair generation "
                      "of heavy parents. Includes arbitrary stable charge ratio, final "
                      "Pauli neglected, two net counted decays per pair. Not a source "
                      "history or a bound for other number-changing channels.")


if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
