"""Primary CLEO III measured continuum Born R_had in 6.964..10.538 GeV.

D. Besson et al. (CLEO), Phys. Rev. D 76, 072008 (2007),
doi:10.1103/PhysRevD.76.072008
Table VII of CLNS-07/1998:
https://www.lns.cornell.edu/public/CLNS/2007/CLNS07-1998/clns-r-sep04.pdf

R(s) = radiatively corrected sigma0(e+e−->hadrons) /
                    (4 pi alpha²/(3s)).
Uncertainties quoted separately as statistical, correlated common
systematic across energies, uncorrelated systematic. A 9.996-GeV
measurement is NOT an exact 10.000-GeV measurement; NEVER silently
replace continuous R(s), interpolate across nearby Upsilon resonances
or claim invisible vector BR known without the full dark+SM channels.
The table has precisely the seven published energy points.
"""
from __future__ import annotations
import math


PRIMARY_DOI="10.1103/PhysRevD.76.072008"
PRIMARY_URL="https://doi.org/"+PRIMARY_DOI
PAPER_TABLE="Table VII, CLEO III continuum Born cross sections"
# Energy; R; stat; common systematic; uncorrelated systematic.
CLEO2007_DATA=(
    (10.538,3.591,.003,.067,.049),
    (10.330,3.491,.006,.058,.055),
    (9.996,3.497,.004,.064,.043),
    (9.432,3.510,.005,.066,.037),
    (8.380,3.576,.024,.058,.025),
    (7.380,3.550,.019,.058,.020),
    (6.964,3.597,.033,.057,.020),
)


def measured_R_at_or_near(mV=10.,max_gap_GeV=.01):
    """Nearest PUBLISHED CLEO continuum R: report mismatch, do not interpolate.

    This is measured continuum R at an energy near mV; extrapolation
    within max_gap is a model ASSUMPTION. At very narrow mV width the
    4MeV gap can matter. R is not an in-medium vector spectral function.
    """
    if not (math.isfinite(mV) and mV>0 and max_gap_GeV>=0):
        raise ValueError("Positive finite mass and nonnegative energy gap")
    x=min(CLEO2007_DATA,key=lambda z:abs(z[0]-mV))
    delta=abs(x[0]-mV)
    if delta>max_gap_GeV:
        raise ValueError("No CLEO 2007 point within stated mass gap")
    E,R,stat,common,uncorr=x
    sigma=math.sqrt(stat*stat+common*common+uncorr*uncorr)
    return dict(experiment="CLEO III",
                physical_energy_GeV=E,target_vector_mass_GeV=mV,
                energy_difference_GeV=mV-E,
                mass_gap_GeV=delta,
                R_had_continuum_Born=R,
                stat_error=stat,
                correlated_systematic_error=common,
                uncorrelated_systematic_error=uncorr,
                total_one_point_quadrature_error=sigma,
                relative_one_point_error=sigma/R,
                original_units="R is ratio of radiatively corrected continuum hadron Born cross section to massless muon point cross section",
                DOI=PRIMARY_DOI,
                source=PRIMARY_URL,
                table=PAPER_TABLE,
                experimental_R_at_nearby_not_exact_mass=delta!=0,
                continuum_threshold_resonance_and_vacuum_only=True,
                no_dark_vector_model_fit=True)


if __name__=="__main__":
    import json
    print(json.dumps(measured_R_at_or_near(),indent=2))
