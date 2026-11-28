# Tabulated late radiation and a check of the residual vacuum source

Prepare the [external central EOS files](INDEPENDENT-BENCHMARK-SETUP.md).
The wider table supplies photon-temperature-dependent effective entropy
and radiation density through the charged lepton transitions. Interpolate
the positive g_s and g_r in log(T), solve s(T)b³=s(T_initial), and use
rho_rad=(pi²/30)g_r(T)T⁴. This is an effective adiabatic radiation history,
not a separate neutrino transport or visible portal backreaction solver.
It avoids using the artificially slow freely redshifted photon-only
comparison as a cosmological prediction.

With the starting temperature near 12 MeV, at b=410 the photon
temperature is 40.94 keV, versus 29.27 keV for freely redshifted photons.
The radiation-only H is 2.5465 times the earlier photon-only reference.
A kinetic run can optionally include its retained dark density on top
of this radiation; neither this table nor the simple comparator
describes unspecified exterior dark matter.

The independent prescribed-background vacuum production integral between
12 and 1 MeV adds about 0.772% to the direct light yield and 3.62e-6 of
the heavy-prompt yield computed at 12 MeV, relative to each prior total.
Within 4 to 1 MeV the increment is 8.6e-6 of the direct yield.
This quantifies omitted late forward production for that comparator.
It does not bound inverse portal reactions, model how newly injected
quanta feed the lowest-momentum bins, or establish the full relic yield.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1
python -W error late_radiation_background_reference.py EOS.late.dat
python -W error residual_visible_pair_production_reference.py EOS.late.dat results/visible-production-cumulative-2026-11-26.json
```
