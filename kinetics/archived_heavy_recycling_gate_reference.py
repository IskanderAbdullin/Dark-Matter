"""Pauli-free heavy-parent production on represented fixed input supports."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import cropped, gross, m, M, mp, points, snapshot_functions, y
from late_radiation_background_reference import TabulatedRadiation, MPL


def angular_integral(s):
    p=np.sqrt(np.maximum(0.,s/4-m*m)); k=np.sqrt(np.maximum(0.,s/4-M*M))
    span=2*p*k
    c=s/2-M*M-m*m
    active=s>=4*M*M
    # Rational endpoint avoids subtracting c and span at large energy.
    zp=np.where(active,mp*mp+(M*M-m*m)**2/np.where(active,c+span,1.),1.)
    zm=zp+2*span
    B=(M+m)**2-mp**2
    ratio=np.divide(np.log1p(2*span/zp),span,out=np.asarray(2/zp).copy(),where=span>0.)
    I=2*(1+B*B/(zp*zm))+2*B*ratio
    return np.where(active,I,0.)


def direct_angular_check():
    x,w=roots_legendre(96)
    err=[]
    for s in np.geomspace(4*M*M*1.00001,.9,30):
        c=s/2-M*M-m*m
        span=2*np.sqrt(s/4-m*m)*np.sqrt(s/4-M*M)
        lo,hi=np.log([mp*mp+c-span,mp*mp+c+span])
        z=np.exp((hi+lo)/2+(hi-lo)/2*x)
        transfer=mp*mp-z
        val=(hi-lo)/(2*span)*np.sum(w*z*(((M+m)**2-transfer)/z)**2)
        err.append(abs(angular_integral(np.array(s))/val-1.))
    assert max(err)<2e-9
    return max(err)


def vacuum_upfold(edges,occ,nx,nz):
    p,e,w=points(edges,m,occ,nx)
    E=e[:,None]; F=e[None,:]
    P=p[:,None]*p[None,:]
    cap=np.clip((E*F+m*m-2*M*M)/P,-1.,1.)
    x,z=roots_legendre(nz)
    mu=((cap-1)/2)[:,:,None]+((cap+1)/2)[:,:,None]*x
    s=2*m*m+2*(E[:,:,None]*F[:,:,None]-P[:,:,None]*mu)
    k=np.sqrt(np.maximum(0.,s/4-M*M))
    rel=y**4*2*k/np.sqrt(s)/(64*np.pi*(E*F)[:,:,None])*angular_integral(s)
    weight=(w[:,None]*w[None,:]*(cap+1)/4)[:,:,None]*z
    # For distinct opposite-charge pairs, each leg has g=2.
    return float(4*np.sum(rel*weight))


def run(table,paths,crop_last):
    cases=[(Path(p).name,json.loads(Path(p).read_text())["snapshots"][-1]) for p in paths]
    if crop_last:
        tag,base=cases[-1]
        cases.append((tag+f": same-state crop to q/T_R={crop_last:g}",cropped(base,crop_last)))
    rows=[]
    for name,s in cases:
        p=np.array(s["lab_momentum_edges_GeV"])
        f,_,g,_=snapshot_functions(s)
        T=s["T_visible_GeV"]
        H=np.sqrt(TabulatedRadiation(table,T).state(1.)[1]/3)/MPL
        nt=np.sum(s["represented_species_count_a3"])/s["scale_factor"]**3
        exact=[]
        for order,chi in [((3,8,8,3),(4,12)),((5,12,12,4),(7,21)),((8,16,16,5),(10,28))]:
            _,ev,_=gross(p,mp,M,f,np.zeros_like,True,*order,exchange_mass=m)
            a=2*ev[1]/(nt*H)
            b=2*vacuum_upfold(p,g,*chi)/(nt*H)
            exact.append(dict(momentum_and_angle_order_phi=order,radial_angle_order_chi=chi,
                forward_phi_phi_heavy_count_over_total_NH=a,
                forward_opposite_charge_chi_heavy_count_over_total_NH=b,
                upper_including_both_same_charge_heavy_channels_over_total_NH=a+3*b,
                mean_heavy_pair_event_energy_from_scalars_GeV=ev[2]/ev[1] if ev[1]>0. else None))
        cutoff_phi=p[-1]/np.sqrt(M*M-mp*mp)
        cutoff_chi=p[-1]/np.sqrt(M*M-m*m)
        if cutoff_phi<=1.:
            assert max(q["forward_phi_phi_heavy_count_over_total_NH"] for q in exact)==0.
        if cutoff_chi<=1.:
            assert max(q["forward_opposite_charge_chi_heavy_count_over_total_NH"] for q in exact)==0.
        rows.append(dict(snapshot=name,T_GeV=T,orders=exact,
            conditional_redshift_only_closure_b_phi_phi=cutoff_phi,
            conditional_redshift_only_closure_b_stable=cutoff_chi))
    return rows


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("paths",nargs="+")
    p.add_argument("--crop-last",type=float)
    a=p.parse_args()
    print(json.dumps(dict(independent_opposite_charge_angle_relative_error=direct_angular_check(),
        cases=run(a.table,a.paths,a.crop_last),
        scope="Fixed represented incoming momentum cells, unpolarized and charge-symmetric. Neglect final Pauli to overestimate gross rates. Prompt decay of all newly made heavy parents would create one net dark particle per parent. No new heavy inverse, time-dependent redistribution, stimulation of decay, or unlimited momentum support is inferred."),
        indent=2,allow_nan=False))
