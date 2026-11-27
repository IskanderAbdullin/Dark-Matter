# The cut real-radiative correction does not remove the bare native sign witness

Extend the [three real-photon two-to-two crossings](INDEPENDENT-EARLY-POLE-YUKAWA-REAL-PHOTON-CROSSINGS-2026-11-27.md) by the vacuum real-radiative decay and its inverse, P ↔ b + scalar + γ. The photon now leaves the decay. The two same-charge diagrams are crossed coherently. A signed-photon matrix trace in a separately constructed parent-rest-frame three-body phase space agrees with the reduced invariant expression, and the associated quantum emission and inverse weights satisfy the common-temperature equilibrium gate. In rate integration the parent mass shell fixes the photon direction relative to the plasma-frame daughter–scalar pair; all allowed energies in the fixed photon interval [0.5,4] are integrated. Both fermion momenta must remain in [0.45,5]. These cuts define a finite observable, not an infrared-inclusive correction.

At q=4.1e−6 on the unmodified archives, the new stimulated terms are +3.80871 (fine) and +3.63269 (third), omitting the common g²y². The 10→16 quadrature differences are below 0.01%. Summing these with the previously tested three finite real-photon crossings, the stimulated term at g=.005 is 7.96e−6 and 1.014e−5 relative to the corresponding bare direct rate, respectively.

For the bounded parent mixtures from the [native sign witness](INDEPENDENT-EARLY-POLE-NATIVE-APERTURE-WITNESS-2026-11-27.md), keep the same scalar and daughter distributions, and evaluate all four photon contributions on the exact native substeps. At this single fine mode the bare direct rate at 85% and 97% mixtures is −11.4494 and −24.5112 (omitting y²); the added real-photon terms in those same units are +.00067645 and +.00068084. The linear mixture fraction at which this bare mode crosses zero moves from .74481384 to .74482002, a change of about 6.18e−6 in mixture fraction. At this same q the third bare rate for these witnesses is positive; this is compatible with the opposite sign in the earlier whole-grid witness because the sign correction depends on the scalar node. No sign changes in the finite real correction at the tested points.

No charged or photon dispersion corrections, external soft photons, virtual amplitudes, coherent multi-scattering, or evolution are included. In particular the cut real correction is not a complete result at order g²y² and is not a bound on the dressed source near a cancellation.

From `kinetics/` with `OPENBLAS_NUM_THREADS=1`:

```sh
for g in fine third; do
  for n in 5 10 16; do
    python early_pole_frozen_yukawa_radiative_decay_reference.py "$g" --order "$n"
  done
  for w in .85 .97; do
    python early_pole_frozen_yukawa_photoproduction_source_reference.py "$g" --order 10 --witness-mix "$w"
    python early_pole_frozen_yukawa_photon_pair_crossing_reference.py "$g" --order 10 --witness-mix "$w"
    python early_pole_frozen_yukawa_radiative_decay_reference.py "$g" --order 10 --witness-mix "$w"
  done
done
```

Raw: `results/early-pole-frozen-yukawa-real-radiative-and-witness-2026-11-27.json.gz`.
