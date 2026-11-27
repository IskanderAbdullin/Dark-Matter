"""Kinetic marks of the actually evolved retained cold distributions."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_angular_mode_initializer import regrid
from early_pole_frozen_dynamic_vlasov_reference import load


def build(archive,angular=None):
    with gzip.open(archive,'rt') as f:
        v=json.load(f)['branches']['steps']
    s=v['snapshots'][0]
    hard=v.get('aperture_subdivision',1)
    C=Cold(v['grid'],1,0,len(s['cold_bins']),True,exposure=v['duration'],
           reference_aperture_grid=v['aperture_grid'],archive_projection='steps',
           aperture_subdivision=hard)
    np.testing.assert_array_equal(C.M.init,s['hard_moments'])
    np.testing.assert_array_equal(C.ns,s['cold_bins'])
    if angular is not None:
        r=load(f"results/early-pole-spectrum-{v['grid']}-m06-2026-11-26.json.gz")
        C=regrid(C,r,angular['shared_scalar_volume_edges'],2)
    return C,hard


def diagnostic(C,hard,raw):
    a,b=raw
    mode=[]
    for x,label in zip(raw,('native','local-daughter'),strict=True):
        if (x.get('prepared_charged_swap','native')!=label or x.get('aperture_subdivision',1)!=hard
                or any(x[k]!=a[k] for k in ('duration','cold_momentum_nodes_over_T','history_exposure','rtol','base_atol'))
                or not x.get('elastic_channels',True) or x.get('selected_elastic_channel','both')!='both'
                or abs(x['energy_defect'])>=3e-12 or abs(x['entropy_defect'])>=3e-11
                or abs(x['hard_number_defect'])>=3e-12):
            raise ValueError('Controls or original solver gates differ.')
        z=np.asarray(x['cold_mode_count_history'])
        np.testing.assert_array_equal(x['cold_momentum_nodes_over_T'],C.q)
        np.testing.assert_array_equal(z[0],C.ns)
        np.testing.assert_array_equal(z.T.copy(order='C').sum(axis=0),x['cold_count_history'])
        if z.shape!=(len(x['history_exposure']),len(C.q)) or not np.all(z>0.):
            raise ValueError('Invalid number history.')
        mode.append(z)
    x,y=mode
    p=C.q
    masses=[C.scalar(z)[-1] for z in (x[0],x[-1],y[-1])]
    mi=0.5*(masses[1]+masses[2])
    A,B=[np.hypot(p,m) for m in masses[1:]]
    weights=p*p/(3*A)
    Aex,Bex=[np.hypot(C.spq,m) for m in masses[1:]]
    bracket=1.+.5*C.lam*(float(y[-1]@(1./(A*B*(A+B))))
                              +float(C.spect@(1./(Aex*Bex*(Aex+Bex)))))
    pole_sq_difference=.5*C.lam*float((y[-1]-x[-1])@(1./A))/bracket
    inferred_gap_difference=pole_sq_difference/(masses[1]+masses[2])
    if abs(inferred_gap_difference-(masses[2]-masses[1]))>5e-17:
        raise ValueError('Independent subtracted gap equation disagrees with direct roots.')
    pole_weights=-p*p*pole_sq_difference/(3*A*B*(A+B))
    # Difference of q^2/(3E) formed without subtracting two nearly equal weights.
    groups=np.searchsorted(np.geomspace(1e-6,1e-5,9)**3,p**3)-1
    off={}
    for k,N,L in (('increment',x[-1]-x[0],y[-1]-y[0]),('terminal',x[-1],y[-1])):
        d=y[-1]-x[-1]
        if not np.all(np.isfinite(N)&(N>=0.)) or not np.all(np.isfinite(L)&(L>=0.)):
            raise ValueError('This diagnostic requires nonnegative net increments in both arms.')
        NN=fsum(N); DD=fsum(d)
        reference=float(weights@N)
        W=reference/NN
        n_term=DD*W
        s_term=float((weights-W)@d)
        m_term=float(pole_weights@L)
        direct=float((weights+pole_weights)@L-reference)
        total=n_term+s_term+m_term
        if abs(total-direct)>16*np.spacing(reference):
            raise ValueError('Stable pressure decomposition did not reproduce the direct difference.')
        pg=np.bincount(groups,weights=weights*d,minlength=8)
        ng=np.bincount(groups,weights=d,minlength=8)
        Q=float(p@N)/NN
        s2=float((p-Q)@d)/(NN+DD)
        lk=float((p*p/(np.hypot(p,mi)+mi))@N)
        off[k]=dict(native_count=NN,paired_count=DD,native_kinetic_pressure=reference,
                    stable_paired_pressure=total,relative_paired_pressure=total/reference,
                    count_term=n_term,shape_term=s_term,pole_term=m_term,
                    count_term_rel=n_term/reference,shape_term_rel=s_term/reference,pole_term_rel=m_term/reference,
                    paired_group_count=ng.tolist(),paired_group_pressure_at_native_pole=pg.tolist(),
                    paired_pressure_nonzero_groups=[int(i+1) for i in np.flatnonzero(pg!=0.)],
                    paired_conditional_mean_p_shift=s2,native_retained_kinetic_energy_at_midpoint_pole=lk)
    scalar=[]
    for count,m in zip((x[-1],y[-1]),masses[1:],strict=True):
        wave=np.r_[p,C.spq]
        num=np.r_[count,C.spect]
        def state(z):
            q=wave*z
            c=num*z**3
            gap=lambda u:u*u-.5*C.lam*(c@np.reciprocal(np.hypot(q,u)))
            mm=brentq(gap,1e-8,.03,xtol=1e-17)
            E=np.hypot(q,mm)
            V=mm**4/(2*C.lam)
            return mm,float(c@E-V),float(c@(q*q/(3*E))+V),V
        m0,rho,P,V=state(1.)
        if abs(m0-m)>3e-17 or abs(rho-C.scalar(count)[1])>3e-17:
            raise ValueError('Independent conformal gap did not match the archived chart.')
        h=1e-3
        rr=[state(np.exp(k*h))[1] for k in (-2.,-1.,1.,2.)]
        Pd=(rr[0]-8*rr[1]+8*rr[2]-rr[3])/(36*h)-rho
        if abs(Pd-P)>2e-9*abs(P):
            raise ValueError('Independent uniform-volume pressure derivative disagrees.')
        scalar.append(dict(energy=rho,pressure=P,pressure_to_energy=P/rho,
                     relative_trace=(rho-3*P)/rho,pressure_from_dilation=Pd,
                     hartree_pressure_fraction=V/P))
    vc=pole_sq_difference*(masses[1]**2+masses[2]**2)/(2*C.lam)
    exterior_kinetic_delta=-pole_sq_difference*float(C.spect@(C.spq**2/(3*Aex*Bex*(Aex+Bex))))
    deltaP=off['terminal']['stable_paired_pressure']+exterior_kinetic_delta+vc
    deltaR=float((y[-1]-x[-1])@A)+pole_sq_difference*(
            float(y[-1]@(1./(A+B)))+float(C.spect@(1./(Aex+Bex))))-vc
    if abs(deltaR-3*deltaP)>1e-16*scalar[0]['energy']:
        raise ValueError('Subtracted full conformal scalar stress does not close to trace zero.')
    return dict(scope='Evaluated scalar-only stress in the fixed held conformal Hartree model; no gravity or inferred physical relic.',
                conformal_scalar_arms=scalar,paired_conformal_scalar=dict(
                    retained_kinetic_pressure_term=off['terminal']['stable_paired_pressure'],
                    exterior_kinetic_pressure_term=exterior_kinetic_delta,
                    hartree_pressure_term=vc,total_pressure=deltaP,total_energy=deltaR,
                    relative_pressure=deltaP/scalar[0]['pressure'],
                    relative_energy=deltaR/scalar[0]['energy']),
                hard_depth=hard,cold_modes=len(p),initial_gap=masses[0],native_final_gap=masses[1],
                local_final_gap=masses[2],
                stable_gap_difference_from_subtracted_equations=inferred_gap_difference,
                independent_gap_difference_discrepancy=float(inferred_gap_difference-(masses[2]-masses[1])),
                native_and_local=off)


def r(path):
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as f: return json.load(f)
    return json.loads(Path(path).read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','native','local'):p.add_argument(k)
    p.add_argument('--angular-mesh')
    a=p.parse_args()
    C,h=build(a.archive,r(a.angular_mesh) if a.angular_mesh else None)
    print(json.dumps(diagnostic(C,h,[r(a.native),r(a.local)]),indent=2,allow_nan=False))
