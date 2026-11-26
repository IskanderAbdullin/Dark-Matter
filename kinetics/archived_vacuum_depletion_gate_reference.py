"""Optimistic vacuum-only depletion after actual retained-count checkpoints."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import C, MPL
from three_species_decay_onsager_reference import M, m, mp, y


def comparison(radiation, normalization, shapes, target=.120):
    fac=conversion(radiation)["Omega_h2_per_GeV_times_yield"]
    table=TabulatedRadiation(radiation,.05)
    t,gs,_=np.loadtxt(radiation).T
    # Upper across 0 <= T <= 50 MeV, including the T=0 endpoint. PCHIP
    # of log(g_s) has no overshoot between knots. Extend by constant g_s
    # below the first positive table knot; its value equals the T=0 row.
    assert t[0]==0. and gs[0]==gs[1]
    x,w=roots_legendre(12)
    K=y**4/(8*np.pi*M*M)
    a_gamma=np.sqrt(2*C/3)/MPL
    out=[]
    for row,path in zip(normalization,shapes,strict=True):
        snap=json.loads(Path(path).read_text())["snapshots"][-1]
        end=row["final_visible_temperature_GeV"]
        assert abs(snap["T_visible_GeV"]/end-1.)<1e-8
        gmax=float(gs[:np.searchsorted(t,end,side="right")+1].max())
        # Total stable particle+antiparticle Y, for a charge-symmetric gas.
        Y=row["yield_per_retained_species"][1]
        initial=fac*m*Y
        assert abs(initial/row["conditional_Omega_h2"][1]-1.)<2e-14
        optical_upper=K*Y/2*(4*C/3)*gmax*end/a_gamma
        # Sharper control uses a radiation-only H, still deliberately
        # dropping matter and all post-checkpoint additional injection.
        upper=np.log(end)
        lower=table.lo
        edges=np.unique(np.r_[lower,np.clip(table.sg.x,lower,upper),upper])
        integ=0.
        for l,r in zip(edges,edges[1:]):
            u=(l+r)/2+(r-l)/2*x
            H=np.sqrt(C/3)*np.exp(2*u+table.rg(u)/2)/MPL
            s=4*C/3*np.exp(3*u+table.sg(u))
            integ+=(r-l)/2*np.sum(w*s/H*(3+table.sg.derivative()(u))/3)
        g_r=float(np.exp(table.rg(lower)))
        integ+=(4*C/3)*gs[1]*np.exp(lower)/(np.sqrt(C*g_r/3)/MPL)
        optical_rad=K*Y/2*integ
        assert 0.<optical_rad<optical_upper
        pmax=snap["lab_momentum_edges_GeV"][-1]
        # Each final scalar lies above this energy on the initial retained
        # incoming support; one of the two necessarily has E >= m.
        max_Etot=2*np.hypot(m,pmax)
        Emin=m*m/max_Etot
        out.append(dict(final_visible_temperature_GeV=end, stable_fixed_count_Omega_h2=initial,
            upper_optical_depth_photon_only=optical_upper, radiation_table_optical_depth=optical_rad,
            lower_Omega_h2_photon_bound=initial/(1+optical_upper),
            lower_Omega_h2_radiation_table=initial/(1+optical_rad),
            necessary_final_F_if_stimulation_is_uniformly_bounded_photon_only=max(0.,np.sqrt((initial/target-1.)/optical_upper)-1.),
            necessary_final_F_if_stimulation_is_uniformly_bounded_radiation_table=max(0.,np.sqrt((initial/target-1.)/optical_rad)-1.),
            retained_fermion_max_lab_momentum_GeV=pmax,
            conservative_scalar_min_lab_energy_GeV=Emin,
            corresponding_scalar_min_lab_momentum_GeV=np.sqrt((Emin-mp)*(Emin+mp))))
    assert min(r["lower_Omega_h2_photon_bound"] for r in out)>target
    return dict(vacuum_sigma_v_upper_GeV_minus2=K, benchmark_target=target, rows=out,
        scope="Symmetric retained stable archive counts, adiabatic visible entropy and expanding radiation floor; only vacuum chi2-antichi2 to two phi, with deliberately maximal removal and no source after the archive. Neither general model exclusion nor a prediction with stimulation or other loss/dilution.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("radiation")
    p.add_argument("normalization")
    p.add_argument("shapes", nargs="+")
    a=p.parse_args()
    rows=json.loads(Path(a.normalization).read_text())["at_fixed_archive_counts"]
    print(json.dumps(comparison(a.radiation,rows,a.shapes),indent=2,allow_nan=False))
