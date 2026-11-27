# Fold the measured charged-pion window into the same cold recoil cut

The [charged-lepton comparison](INDEPENDENT-WARM-LATE-COLD-VISIBLE-RECOIL-ENCLOSURE-2026-11-26.md)
left the hadronic source separate. Add only the exclusive
[KLOE combined vacuum charged-pion band](https://arxiv.org/abs/1711.03085),
using the same 12–50 MeV visible history, p/T<.1 at reference
equality, and unpolarized open decay. Use the tabulated
vacuum cross sections as piecewise constant on their
published bins, .10≤s≤.95 GeV². No extrapolation is used.
The electron/muon source and the incoming charged-pion
source are distinct; the pion input uses the same vacuum
spectral transfer convention as the existing [KLOE fold](kloe_pion_band_portal_reference.py).

For every invariant-mass node solve the exact minimum
total energy for one parent to pass the scalar recoil
gate. Include incoming pion Bose distributions in the
bath and the spin-zero incoming pion current's angular
weight, proportional to 1−beta_parent² cos² theta in the
pair rest frame. A direct azimuth average validates the
Legendre reduction to 4.5e−16 absolute. The independent
uncut Maxwell–Boltzmann boost integral recovers its
Bessel-function normalization to 1.2e−15 relative.

The measured-window pion yield is 4.94437e−10 with Bose
weights, or 11.83% of the same held q64 cold yield. It
is 4.92375e−10 in the MB comparison. The quadrature
controls agree to 6e−6 relative. Summing the individual
published marginal uncertainties linearly changes this
yield by 1.16e−11 in the fixed-bin model; this is a
sensitivity check, not a confidence interval or a
treatment of within-bin form-factor shape uncertainty.

Adding these two disjoint numerical vacuum sources gives
about 29.9% of the held cold comparison, before any
unmeasured bands or other sources. The interval bound
for the lepton source is not a bound on this combined
model. This measured exclusive window omits the pion
threshold region below the table, higher invariant masses,
other channels, rescattering, and finite-temperature
spectral changes. It is not a full hadronic estimate.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_visible_kloe_recoil_reference.py EOS.dat \
  results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz KLOExSec.data
~~~

[Rate, data hash and independent controls](results/warm-late-cold-visible-kloe-recoil-2026-11-26.json).
