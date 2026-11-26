# Instantaneous and finite reheating on exactly the same metric

This extends the [separate CDM-driver survey](INDEPENDENT-FINITE-COOLING-CDM-LOS-SURVEY-2026-11-26.md)
with the instantaneous prompt-plus-first-order pair spectrum. The
finite cooling clock is the same stipulated `T ∝ a^-3/8` example
with `Tmax/TR = 5`, not an inferred thermal history. For the 86 eV
scalar the joined inputs have respectively
`Omega h² = 3.04510e-10` (instantaneous) and `2.74200e-9`
(finite), as recorded by the production reference. The 1 MeV
anchor RMS `p/T` is 7.37380 and 2.04729, respectively. Unlike a
density-matched cosmological comparison, these shaped inputs retain
their own tiny absolute counts. Overall normalization cancels in
each scalar density contrast on the imposed metric.

## Matched probes

Evaluate both shapes with the same mass, metric histories, early
perturbation prescription, momentum/time settings, and exact target
modes. The comparison verifies every native file hash and both
the start and final conformal times through the recorded metadata.
Entries below are the **scalar's** joined density contrast relative
to the driver CDM density contrast, not a total-matter power ratio.

| z | k [h/Mpc] | Instantaneous | Finite | Paired difference |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 22 | 0.291009 | 0.856585 | 0.565576 |
| 5 | 22 | 0.111030 | 0.736954 | 0.625924 |
| 20 | 4 | 0.369240 | 0.881305 | 0.512066 |
| 20 | 10 | 0.137372 | 0.759673 | 0.622301 |
| 20 | 22 | 0.041915 | 0.604525 | 0.562610 |

The finite shape has the higher joined ratio at all 21 matched
samples. Across them the max refinement difference is `3.4e-8`
and the discrepancy from the energy-weighted prompt-plus-recoil
calculation is `7.1e-9`. The instantaneous joined half-amplitude
brackets are `[10,16]` at z=0, `[4,8]` at z=5, and `[2,4]`
at z=20, all in h/Mpc. The z=20 finite case is still above half
through the common grid at k=22; the separately sampled extension
brackets this level in `[34,40]`. These are measured brackets at
sampled points only, not certified continuous thresholds. No
finite high-k extension at z=0 or 5 is used or needed for the
matched comparison.

The paired change shows sensitivity to the production history
in the frozen scalar perturbation **given both source spectra**.
It is not evidence for an abundance fit, a nonlinear halo effect
or a constraint on real matter power. The comparison includes
neither late collisions nor scalar gravitational backreaction,
and does not assign the effect to just one part of the different
early source histories.

## Reproduction

Keep the perturbation files from a single CDM-only native driver.
Generate the instantaneous/finite prompt and joined unremapped
sigma=0.14 tables from the prior
[production spectrum](INDEPENDENT-FINITE-COOLING-CDM-LOS-SURVEY-2026-11-26.md).
Build both recoil tables directly from the detailed pair record:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_matched_history_los_reference.py recoil \
  EOS.dat results/finite-cooling-during-production-pair-shape-2026-11-26.json recoil >recoil.json
~~~

The two rebuilt recoil files reproduce the already used tables
byte for byte. For each respective prompt/joined/recoil triplet and
each z=0,5,20, run the previous LOS command with the **same** raw
metric filenames and report JSON. Finally check the paired invariants:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_matched_history_los_reference.py compare \
  --finite results/finite-cooling-cdm-los-z{0,5,20}-2026-11-26.json \
  --instant results/instant-reheat-cdm-los-z{0,5,20}-2026-11-26.json \
  results/finite-cooling-cdm-los-high-z20-2026-11-26.json \
  results/finite-cooling-matched-history-recoil-2026-11-26.json \
  results/finite-cooling-pair-mixture-export-2026-11-26.json
~~~
