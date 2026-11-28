# A dilute scalar-scattering clock does not justify neglecting scattering

Evaluate a deliberately dilute reference clock on the
[retained early distributions](INDEPENDENT-EARLY-POLE-SPECTRAL-SNAPSHOT-2026-11-26.md).
For a real scalar with amplitude λ and the same momentum-independent
Hartree pole used for the dispersion, assign the identical-particle
Born cross section λ²/(32πs). Average cross section times the Møller
flux factor over the angle between the two incoming momenta and
over both stored shell inventories. The nonrelativistic isotropic
relative-speed formula independently verifies the angular treatment.

At the illustrative parent mass 86 GeV, corresponding to an
86 eV bare scalar for m₀/M=10⁻⁹, the number-weighted dilute event
clock is 0.65–0.79 times radiation-reference H across the third-grid
x=0.05 cases. For daughter/parent mass 0.06 at x=0.10 it is about
1.53–1.58 times reference H on the two extension grids.
At fixed dimensionless trajectory the clock in H units scales as
the inverse of the assumed dimensional parent mass. These are total
dilute event counts per particle, not rates for damping a mode or
transporting momentum in the actual medium.

The same spectra have λ times peak occupancy above one billion
and n/m_eff³ about 2.6 × 10⁷, making the dilute approximation
uncontrolled. The latter scaling is expected because the tadpole
gap approaches m_eff³≈λn/2 in the deep IR. Bose factors,
gain–loss cancellation, medium vertices and correlations could
all change the functional result. This diagnostic gives no
bound on the true transport or condensation timescale, and does
not show that the frozen spectra are observable. An actual
collision or field-evolution treatment is required before neglecting
redistribution or using the conditional equilibrium closure.

~~~sh
python -W error early_pole_scalar_dilute_born_clock_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Raw rates, angle checks and occupancy diagnostics](results/early-pole-scalar-dilute-clock-2026-11-26.json).
