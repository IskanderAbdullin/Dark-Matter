# Chemical number and heat in a finite shared photon reservoir

Promote the [hard central-angle pair check](INDEPENDENT-EARLY-POLE-HARD-GAUGE-PAIR-SOURCE-2026-11-26.md)
to a reversible count-and-heat system for both Dirac species.
This is a fixed-scale thermal closure, not a Friedmann evolution.
Photons are assumed to maintain a Bose shape at zero chemical
potential; the two fermions are assumed to maintain FD shapes
with a common temperature and independent chemical potentials.
Photon thermalization is an assumption, not supplied by the
single pair channel modeled here.

Use exact pair lab energies and finite angular weights from the
vacuum amplitude, with the same hard radial, invariant-mass and
center-of-mass angle cuts. In the lab, the boost-direction cosine
of the incoming photon in the center-of-mass frame is
(p₁−p₂)/|P|. This fixes the azimuthal average for final FD factors.
Reverse/forward for each event is
exp[2μ_f/T_f + E_total(1/T_γ−1/T_f)].
There are two fermion particles per event but only one pair
energy. The photon reservoir loses exactly the energy added
to the fermions in this nonexpanding comparison.

At equal temperature and zero chemical potential, the velocity
and count source vanish pointwise. At unequal temperatures,
one can tune the chemical potential to make the *integrated
count* vanish and still have nonzero heat from the hotter side.
For Tγ=1 and Tf=.8 in reference units, the finer normalized
bare example gives heat 1.693 × 10⁻⁴ at zero count. Reversing
the temperature ordering reverses its sign. The integrated
entropy production is positive. These checks do not rely on
assuming energy deposition proportional to net number.

At fixed early x=.05 and g=.001, the active reservoir case
ends the chosen reference exposure 0.8 with daughter count
.15583 in Tref³ units and Tγ/Tref=.90503 on the finer
quadrature. The thermostated counterpart has count .19499
and needs .37923 Tref⁴ of imposed energy. An initially colder
photon bath instead gains energy as the parent first
annihilates, and reaches daughter count .03701 in this test.
The longer fixed-scale reference approaches the separately
computed zero-chemical-potential temperature .88388.
These exposures are not cosmic times or full in-medium rates.

~~~sh
python -W error early_pole_finite_shared_photon_chemistry_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --order 22 --refined-angle --end .8
python -W error early_pole_finite_shared_photon_chemistry_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --end 12
~~~

[Raw closed, thermostated and zero-count heat checks](results/early-pole-finite-shared-photon-chemistry-2026-11-26.json.gz).
