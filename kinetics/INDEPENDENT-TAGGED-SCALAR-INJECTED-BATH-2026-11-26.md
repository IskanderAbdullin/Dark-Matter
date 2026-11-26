# Test-scalar transport in a frozen, nonthermal injected fermion bath

Here the fermion distribution is taken from the [collisionless injection and prompt-decay spectrum](INDEPENDENT-EMPTY-LEPTON-INJECTION-2026-11-26.md), rather than a thermal fit with the same density or energy. There is no forward feedback on that bath or later cosmological evolution. A test scalar is sufficiently dilute that scalar stimulation is suppressed in its scattering operator. This is a linear single-particle transport diagnostic, not a near-condensate calculation.

Evaluate exact below-pole elastic geometry with nonthermal incoming f_chi(E) and outgoing fermion blocking, instead of assuming reciprocal scalar transitions. A positive finite-volume number generator follows; the bath itself does not satisfy scalar detailed balance, so the ordered matrix is not symmetrized. The cutoff for the auxiliary hitting problem is obtained independently from the **actual specified** f_chi using the momentum-resolved inverse-decay rate. The heavy parent population is set negligible under the prompt-decay premise. The high fermion energies not supported by the injected spectrum are not extrapolated. Independently check mean times using 18,000 continuous-time stochastic paths. Repeat with increased angular/energy orders and with piecewise-constant occupations in place of shape-preserving interpolation.

| Reheat → snapshot visible T (GeV) | First dissipative window starts at p/T | Mean test-particle time to reach it, in Hubble times |
|---|---:|---:|
| .020 → .018 | .287 | 4.07 |
| .050 → .045 | .0395 | .0276 |
| .100 → .090 | .00899 | .00116 |

The mean first-passage time changes by less than 0.75% across the stated profile/degree tests. An absolute direct single-scattering rate from the lowest bin is harder to quadrature-converge: at .10→.09 it moves from .125 to .101 in H units. This does not drive the repeated-path result; do not quote that direct value as converged. At the coldest snapshot, a first-passage time exceeding H^-1 makes a frozen-background lifetime a poor predictor of the actual evolving system. At higher reheating, the static dilute transport can be rapid, but the fermion spectrum may itself change on a similar timescale. Does not measure a condensate, actual population transfer, or an abundance.

Reproduce behind the same central [SM EOS table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat):
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > raw-spectrum.json
python -W error kinetics/tagged_scalar_nonthermal_bath_reference.py EOS_TABLE raw-spectrum.json
```
