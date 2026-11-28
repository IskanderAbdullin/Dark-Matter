"""Separate an unconditional kinetic photon moment from a selected-event spectrum."""

import json
import numpy as np


def fixture(k,T,imbalance,D):
    # N is the initial photon-bin occupancy in a canonical single oscillator.
    # All three channels are positive selected-event intensities, not a
    # perturbative value for a real plasma. The no-transition intensity is C+D*N.
    k,T=np.longdouble(k),np.longdouble(T)
    z=np.exp(-k/T)
    F=np.longdouble('0.7')*z*(1+imbalance)
    R=np.longdouble('0.7')
    C=np.longdouble('2.1')
    D=np.longdouble(D)
    nbar=z/(1-z)
    var=nbar*(1+nbar)
    W=F*(1+nbar)+R*nbar+C+D*nbar
    change=F*(1+nbar)-R*nbar
    cov=(F+R+D)*var
    conditional=(change+cov)/W
    # Truncate a summation only at an exponentially negligible Planck tail.
    N=int(np.ceil(44.*T/k))
    n=np.arange(N+1,dtype=np.longdouble)
    density=(1-z)*np.exp(-k*n/T)
    direct_w=F*(n+1)+R*n+C+D*n
    direct_event=F*(n+1)-R*n
    before=n@density
    old_sel=(n*direct_w)@density/(direct_w@density)
    new_sel=((n*direct_w+direct_event)@density)/(direct_w@density)
    u=np.array([float(direct_event@density),float(new_sel-old_sel),float(new_sel-before)])
    a=np.array([float(change),float(change/W),float(conditional)])
    absolute=np.abs(u-a)
    scaled=absolute/np.maximum(np.abs(a),1.)
    reciprocity=np.max(np.abs((density[:-1]*F*(n[:-1]+1))/(density[1:]*R*n[1:])-1))
    if max(scaled)>5e-10 or abs(reciprocity-abs(imbalance))>6e-14:
        raise ArithmeticError((k,T,scaled,reciprocity))
    return dict(k=float(k),T=float(T),departure_from_adjacent_detailed_balance=imbalance,unshifted_event_slope=float(D),
                mean_photon_count=float(nbar),unconditional_number_clock=float(change),
                mean_per_event_physical_heat=float(k*change/W),
                conditional_post_event_minus_its_selected_initial_energy=float(k*(new_sel-old_sel)),
                conditional_post_event_minus_unselected_thermal_energy=float(k*conditional),
                selection_covariance_energy_term=float(k*cov/W),
                scaled_direct_sum_defect=float(max(scaled)),
                adjacent_detailed_balance_relative_defect=float(reciprocity))


def run():
    rows=[fixture(k,T,s,D) for k,T in ((4e-5,1.),(.000435,1.),(.08,.4),(.6,1.))
          for s,D in ((0.,.19),(-.12,.09),(.17,0.))]
    for b in rows:
        if b['departure_from_adjacent_detailed_balance']==0.:
            if abs(b['mean_per_event_physical_heat'])>1e-14 or b['selection_covariance_energy_term']<=0.:
                raise ArithmeticError('Thermal selection failed the interpretation check')
    return dict(scope='Manufactured one-mode positive event selection; not a finite thermal virtual term or a plasma photon spectrum.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
