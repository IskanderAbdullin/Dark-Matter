# A frozen elastic tangent can miss the short nonlinear directional response

Repeat the [frozen tangent test](INDEPENDENT-EARLY-POLE-SCREENED-FROZEN-TANGENT-2026-11-27.md) for the same prepared 97% hard witness, FD susceptibility metric, chosen perturbation and finite graph. Instead of holding the implicit geometry and base state fixed, propagate positive and negative logit offsets through the actual conservative elastic steps and take a centered difference. The dielectric weights are held in both tests; no slow hard or scalar channel is evolved. The scalar readout is a directional derivative evaluated on the final states, not their production rate.

| Fine macro steps over exposure 1e−8 | Surviving norm, nonlinear | Scalar directional derivative, nonlinear |
| ---: | ---: | ---: |
| 1 | .29643 | −3.953e−4 |
| 2 | .21978 | −2.882e−4 |
| 4 | .17737 | −2.312e−4 |
| 8 | .15756 | −2.050e−4 |
| 16 | .14958 | −1.944e−4 |

The frozen fine comparator gives .65922 and −.0022755. The two centered perturbation amplitudes agree to about 4e−9 absolute in the reported norms, so the fine discrepancy is not explained by that differencing choice. The 8→16 nonlinear changes are still about 5% and are not a time-convergence certificate. On third, the one- and two-step nonlinear norms are .3034 and .2367 versus frozen .2401; the respective source derivatives at 1e−8 are positive, +5.783e−5 and +7.187e−5, compared with frozen +6.214e−5. These short nonstationary controls do not furnish a new spectral damping rate or a plasma-instability conclusion.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
 python -W error early_pole_screened_witness_step_derivative_reference.py "$g"
 python -W error early_pole_screened_witness_step_derivative_reference.py "$g" --macro 2 --exposures 1e-8
done
for n in 4 8 16; do
 python -W error early_pole_screened_witness_step_derivative_reference.py fine --macro "$n" --exposures 1e-8
done
```
Raw: `results/early-pole-screened-witness-finite-step-derivative-2026-11-27.json.gz`.
