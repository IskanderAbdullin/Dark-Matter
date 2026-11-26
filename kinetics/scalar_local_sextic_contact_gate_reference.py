"""Conditional vacuum contribution from a local c6*phi^6/6! contact."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import points,snapshot_functions
from conditional_relic_normalization_reference import conversion
from late_radiation_background_reference import TabulatedRadiation
from scalar_quartic_vacuum_multiplication_bound_reference import global_upper
from sm_dark_entropy_exchange_reference import C,MPL
from three_species_decay_onsager_reference import m,mp

D=576*(4*np.pi)**5
C6=2.12e-16 # Approximate zero-momentum illustrative coefficient from PR1.


def angular_average(E,p,cut):
    ef=E[:,None]*E[None,:]
    pq=p[:,None]*p[None,:]
    d=2*(ef+mp**2)
    h=2*pq
    hi=np.clip((d-16*mp**2)/h,-1.,1.)
    lo=np.full_like(hi,-1.) if cut is None else np.clip((d-cut**2)/h,-1.,1.)
    width=np.maximum(hi-lo,0.)
    mid=(lo+hi)/2
    # Uniform cos(theta), including threshold and optional EFT cutoff.
    return width*((d-h*mid)**2+h*h*width*width/12)/(8*ef)


def test_angle():
    from scipy.special import roots_legendre
    x,w=roots_legendre(16)
    p=np.geomspace(1e-7,.20,20)
    E=np.hypot(p,mp)
    errors=[]
    for cut in (None,.005,.010):
        a=angular_average(E,p,cut)
        for i,k in [(0,1),(1,18),(4,12),(6,6),(8,19),(12,15),(18,18),(19,19)]:
            ef=E[i]*E[k];pq=p[i]*p[k]
            d=2*(ef+mp*mp);h=2*pq
            h1=min(1.,(d-16*mp*mp)/h)
            l1=-1. if cut is None else max(-1.,(d-cut*cut)/h)
            if h1<=l1:
                assert a[i,k]==0.
                continue
            z=(h1+l1)/2+(h1-l1)/2*x
            q=(h1-l1)/2*np.sum(w*(d-h*z)**2)/(8*ef)
            errors.append(abs(a[i,k]/q-1.))
    assert max(errors)<3e-12
    return max(errors)


def check(table,files,cuts):
    rows=[]
    for path in files:
        s=json.loads(Path(path).read_text())["snapshots"][-1]
        T=s["T_visible_GeV"]
        edges=np.array(s["lab_momentum_edges_GeV"])
        fp,_,_,_=snapshot_functions(s)
        H=np.sqrt(TabulatedRadiation(table,T).state(1.)[1]/3)/MPL
        controls=[]
        for n in (6,12):
            p,E,w=points(edges,mp,fp,n)
            energy=w@E
            nsc=w.sum()
            samples=[]
            exact=w@angular_average(E,p,None)@w
            no_cut_ceiling=energy**2+2*mp**2*nsc**2+mp**4*(w@(1/E))**2+(w@(p*p/E))**2/3
            assert 0.<exact<=no_cut_ceiling*(1+1e-12)
            for cut in cuts:
                v=w@angular_average(E,p,cut)@w
                assert 0.<=v<=exact*(1.+1e-12)
                K=C6*C6*v/D
                samples.append(dict(max_cm_energy_GeV=cut,
                    pointwise_phase_space_upper_weight_fraction=v/exact,
                    upper_contact_number_gain_over_H_Nscalar=K/(H*nsc)))
            controls.append(dict(radial_order=n,rows=samples))
        rows.append(dict(snapshot=Path(path).name,T_GeV=T,controls=controls))
    return rows


def needed(table,normalization,cuts):
    t,gs,_=np.loadtxt(table).T
    val=conversion(table)["Omega_h2_per_GeV_times_yield"]
    a=np.sqrt(2*C/3)/MPL
    ratio_lambda=7.65175e-8
    K4=global_upper(ratio_lambda)
    sample=[]
    for r in normalization["at_fixed_archive_counts"]:
        fit=r["final_visible_temperature_GeV"]
        phi,chi=r["conditional_Omega_h2"]
        Y=(phi/mp+chi/m)/val
        J=(4*C/3)*max(gs[:np.searchsorted(t,fit)+1])*fit/a
        requisite=1.-phi/(.99*.120)-chi*mp/m/(.99*.120)
        tau_contact=[]
        for cut in cuts:
            K6=C6*C6*cut*cut/D
            tau=K6*Y*J
            # Below cutoff, amplitudes can interfere constructively.
            # Above cutoff we include only the previous quartic bound;
            # the true nonlocal loop is explicitly not included.
            both=(np.sqrt(K4)+np.sqrt(K6))**2*Y*J
            min_contact=max(0.,np.sqrt(requisite/(Y*J))-np.sqrt(K4))*np.sqrt(D)/cut
            tau_contact.append(dict(max_cm_energy_GeV=cut,contact_growth_optical_upper=tau,
                    combined_with_halo_reference_quartic_upper=both,
                    optimistic_min_abs_c6_GeV_minus2_with_halo_reference_quartic_for_99_percent_mass=min_contact,
                    ratio_of_min_contact_to_illustrative_matching=min_contact/C6))
        sample.append(dict(snapshot=r["archive"],T_GeV=fit,rows=tau_contact))
    return sample


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("normalization");p.add_argument("snapshots",nargs="+")
    a=p.parse_args()
    cut=(.005,.010)
    n=json.loads(Path(a.normalization).read_text())
    print(json.dumps(dict(illustrative_local_c6_GeV_minus2=C6,independent_angle_relative_error=test_angle(),
         represented_rate_controls=check(a.table,a.snapshots,cut),global_optimistic_contacts=needed(a.table,n,cut),
         scope="A local sextic EFT toy coefficient only, below the stated CM cutoff. The true momentum-dependent fermion loop, hard states above cutoff, coherent fields, arbitrary Bose stimulation, additional particles and interactions are not constrained."),
         indent=2,allow_nan=False))
