# Dependence on the width of a continuous two-moment ramp

Change the [exact continuous-ramp](INDEPENDENT-EARLY-POLE-SCALAR-CONTINUOUS-POLE-RAMP-2026-11-27.md) half-width from 3% to 12% of the smaller adjacent native spacing in p³. Keep the same tagged momentum, the same two native cells and their total number and energy, and the same subtraction. No finite-subcell approximation is used. The widths are alternative occupation reconstructions, not changes to the parent propagator.

| Archive, exterior net / Href | R=2e−4, ramp 3% / 12% | R=2e−3, ramp 3% / 12% |
| --- | ---: | ---: |
| Fine | −1.767e−17 / −4.408e−18 | −1.198e−16 / −4.408e−17 |
| Third | +1.767e−19 / +1.050e−17 | +9.824e−18 / +9.998e−17 |

Values are the signed local remainder in the zero-width numerical control. For either archive the sign survives changing the ramp width in these panels, but on third the narrow-window magnitude increases about 59 times. Increasing the invariant and core quadrature orders changes the wider-ramp results by less than 1.3e−5 relative. The wider fine ramp has no sampled incoming topology threshold inside the outer window and its remainder scales nearly linearly with R. The third wide profile has a threshold at |x|=.0008395. For resolved small offsets the continuous squared-pole G/x² is approximately −2.207e−14 and +5.242e−14 respectively. This is sensitivity to the arbitrary shape within moment-constrained cells. It is neither a physical width nor a change to a predicted production rate.

From `kinetics/`:
```sh
for g in fine third; do
 OPENBLAS_NUM_THREADS=1 python -W error early_pole_archive_scalar_pole_continuous_ramp_reference.py "$g" --fraction .12
done
```
Raw: `results/early-pole-archive-scalar-wider-pole-ramp-2026-11-27.json.gz`.
