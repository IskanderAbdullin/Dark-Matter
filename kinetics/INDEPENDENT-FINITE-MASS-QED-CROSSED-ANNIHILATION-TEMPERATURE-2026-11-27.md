# The crossed contribution has a photon-temperature zero with other factors frozen

Use the [fixed-mode crossed comparator](INDEPENDENT-FINITE-MASS-QED-CROSSED-ANNIHILATION-SOURCE-2026-11-27.md)
with both charged histograms, masses, and scalar occupation held fixed.
The anti-daughter occupation is still the explicit stored daughter proxy.
Vary **only** the Bose temperature of the photon bath.

Let $A=f_P f_{\bar b}$, $B=(1-f_P)(1-f_{\bar b})$ and $n_s$ be the
held scalar occupation. Under the positive spin and phase-space
measure, the source divided by $n_s$ is the integral of
$(1+1/n_s)A + n_\gamma[(1+1/n_s)A-B]$.
Since $(1+1/n_s)A-B=f_P+f_{\bar b}-1+f_P f_{\bar b}/n_s$,
a strict upper bound using the maximum stored occupations is
negative in both archives. The photon occupation increases
strictly with temperature. The vacuum gain is positive and the
warm endpoint of the bracket is negative, so the isolated
fixed-background contribution crosses zero exactly once.
Detailed balance of the complete system is **not** implied at this zero.

Store the positive gain and loss integrals at each photon node
before applying the Bose factors. Reweight the same cells to
avoid confusing temperature dependence with support changes.
Deliberately split at the reported photon-energy cuts. Verify the
analytic temperature derivative independently against a centered
temperature difference and compare two radial/angular meshes.

| Archive | Photon temperature at zero | Vacuum net / scalar occupation | Net / scalar occupation at Tγ=1 | Fraction of inverse loss below photon energy 5 at Tγ=1 |
|---|---:|---:|---:|---:|
| fine | 0.49937298 | +0.228621 | -3.575604 | 87.82% |
| third | 0.50481099 | +0.236702 | -3.528597 | 87.76% |

The two tested meshes move the zero by at most 3.57e-11 in this frozen fixture; the difference between archive shapes is much larger.

Related matched thermal sources need charged and scalar populations
to be updated together; this frozen experiment cannot infer the
physical plasma temperature or identify a cosmological epoch.

Replay from kinetics with the pinned spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_crossed_annihilation_temperature_reference.py
~~~
[Receipt](results/finite-mass-qed-crossed-annihilation-temperature-2026-11-27.json).
