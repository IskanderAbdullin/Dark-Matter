# Rare warm scalar nodes dominate higher moments of the saved numerical spectra

Use the [three time-matched scalar endpoint histograms](INDEPENDENT-EARLY-POLE-FINE-NATIVE-SCALAR-STARTUP-CLOCK-2026-11-28.md), all at x=.05 on the same fine two-subcell charge chart. Hold the scalar mass at the S2 gap value .0008731833 and use a piecewise-constant occupation in each actual cell. Partition the retained momentum range at q=1e−5, .0003, .01 and .1. Tables report normalized shares of moments of the count measure `q²f(q)dq/(2π²)`. In addition to count, q and q², evaluate kinetic energy `q²/(sqrt(m²+q²)+m)` and the pressure numerator `q²/sqrt(m²+q²)`; the isotropic factor of 1/3 cancels in reported pressure shares.

| Scalar chart | Count at q ≤ 1e−5 | Count at q ≥ .01 | q² moment at q ≥ .01 | Kinetic energy at q ≥ .01 | Pressure numerator at q ≥ .01 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Unrefined 67 bins | 99.998928% | .000930% | 99.999926% | 99.9160% | 99.8348% |
| S2, 93 bins | 99.998930% | .000929% | 99.999933% | 99.9251% | 99.8528% |
| S4, 145 bins | 99.998930% | .000929% | 99.999934% | 99.9263% | 99.8552% |

In this benchmark a high-momentum fraction of fewer than ten particles per million in number controls nearly all kinetic energy and pressure, while the histogram differences in the cold mode can remain large. The pressure here is an instantaneous functional at the stated mass, not a cosmological propagation result. In particular, cold number dominance and agreement of a few total moments cannot establish the stability of the full distribution or certify it for CLASS. The warm fractions are from the declared retained numerical support; the calculation neither validates an actual high-energy population nor extrapolates any tails.

The first three columns are independently checked against exact polynomial integrals with factored power differences for nearly coincident bin cuts. All five measured columns agree under 24→48 quadrature to better than 5e−16, and the bands recover the original counts. From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_scalar_spectral_moment_budget.py \
  results/early-pole-fine-charge2-native-scalar-start4096-x05-state-2026-11-28.json.gz \
  results/early-pole-fine-charge2-scalar{2,4}-x05-state-2026-11-28.json.gz
```
The [receipt](results/early-pole-endpoint-scalar-spectral-moment-budget-2026-11-28.json) includes exact thresholds, absolute moments and each band fraction.
