# Restricted direct scalar abundance over an instantaneous-reheating ladder

The adopted visible-temperature history is a prescribed adiabatic
radiation clock, started at separate reheating temperatures from
25 to 100 MeV. This comparison asks whether changing that temperature
can close the [absolute direct-production gap](INDEPENDENT-WARM-DIRECT-ABSOLUTE-NORMALIZATION-2026-11-26.md)
in the tested vacuum channels alone. Each produced heavy parent is
credited with one scalar daughter and all such daughters are retained.
This is not a lower bound on its final abundance after absorption,
and it is not a bound on a model with nonlinear recycling or extra
sources. The [existing 55.5–56 MeV comparison](INDEPENDENT-WARM-FROM-EMPTY-TEMPERATURE-SENSITIVITY-2026-11-26.md)
answers a different question: its deliberately loose **upper** also
includes nonlinear recycling in another selected network.

For each temperature panel integrate the positive increment in ln a
exactly, enclose s_visible H by its monotone endpoints, and bound the
remaining Bessel kernel and the invariant prefactors separately over
the bin. The upper uses the existing MB charged-lepton substitution,
pion Bose multiplier, measured window with a stated nuisance choice,
two low-mass positive moment caps, and a conditional high-band R_pi≤1.
For the lower, if E1+E2≥sqrt(s) and both Ei≥m_l, the FD denominator
relative to MB is at most
(1+exp(−m_l/T))(1+exp(−(sqrt(s)−m_l)/T)).
Indeed the sum of the two exponentials is convex in the energy
split at fixed sum and largest at an endpoint; this upper then
decreases as the total energy rises. Charged pion BE is at least MB.
Use max(0, bin−2 quoted marginal errors) in the fixed measured
window and discard unmeasured pion bands and all tails in the lower.
The two fixed-bin perturbations are not a statistical confidence
interval. The three invariant/temperature resolutions all tighten
the numerical endpoints.

A separate lab-energy and angular integration retains full FD and
BE factors, uses central measured pion bins only, and folds down
to 6 MeV (the two enclosures continue to .1 MeV). It supplies a
numerical control, not a replacement for either positive enclosure.

| Reheating (MeV) | Conditional lower | Numerical central | Conditional upper | Central fraction from T > 0.9 T_R |
|---:|---:|---:|---:|---:|
| 25 | 1.6826e-06 | 2.7008e-06 | 3.468e-06 | 58.4% |
| 40 | 3.5735e-05 | 5.4119e-05 | 7.1426e-05 | 43.3% |
| 45 | 6.4945e-05 | 9.7239e-05 | 0.00012848 | 40.7% |
| 50 | 0.00010735 | 0.00015898 | 0.00021004 | 38.8% |
| 55.5 | 0.00017333 | 0.00025345 | 0.0003342 | 37.6% |
| 60 | 0.0002462 | 0.00035603 | 0.00046883 | 36.9% |
| 75 | 0.00065994 | 0.00091985 | 0.0012227 | 35.9% |
| 100 | 0.0021701 | 0.0029167 | 0.0041955 | 33.5% |

The 50 MeV numerical central total reproduces the earlier result to
about 8e−7 relative, despite integrating through the colder tail. Across
this ladder the coarser/finer angular and time fold disagrees by no
more than 4.3e−6 relative per included species. In the 50 MeV case the
independent enclosure certifies at least 20.9% of the yield above
45 MeV under all the stated channel and positive cap conditions,
compared with 38.8% in the central quadrature. This is an integral
residency statement, not production at a single temperature.

Even at 100 MeV the direct yield upper is 0.00420, less than 7%
of .060 (half of the assumed total .120). Neither this result nor
the numerical central value constrains dark recycling, an arbitrary
thermal hadronic spectrum, the missing stable species, nonstandard
expansion, or reheating outside the tested range.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_temperature_bracket_reference.py EOS.late.dat KLOExSec.data \
  "$D/warm-late-cold-pion-dispersive-recoil-2026-11-26.json" \
  "$D/warm-direct-vacuum-absolute-source-envelope-2026-11-26.json"
~~~

The [machine-readable ladder](results/warm-direct-vacuum-temperature-bracket-2026-11-26.json)
pins both source files by SHA-256. EOS and KLOE data preparation are
described in the [benchmark setup](INDEPENDENT-BENCHMARK-SETUP.md)
and the [measured-band calculation](INDEPENDENT-WARM-LATE-COLD-VISIBLE-KLOE-RECOIL-2026-11-26.md).

The [temperature-dependent spectra](INDEPENDENT-WARM-DIRECT-TEMPERATURE-SHAPES-2026-11-26.md)
retain source attribution and test the redshifted momenta independently of this count bracket.
