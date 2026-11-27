# Preserve the aperture-local input while matching total daughter moments

Keep the two continuous daughter cells inside the initial cold aperture
exact, while applying two common affine FD logit shifts to the other
18 daughter cells. The two free coefficients match the native total
retained daughter number and bare energy. Other hard species, the cold
seed and all support endpoints are unchanged. This manufactured match
does not hold the mean-field energy or all screening moments fixed.

The raw local substitution changes retained daughter number by
-8.81e-5 and bare energy by -1.29e-4. Retuning leaves defects
(-3.47e-18, 0), and the initial mean-source shift changes from
-0.01551056101 to -0.01551107287. Its magnitude differs by about
33 ppm from the unmatched local shift. This is only a static source
control. The reverse experiment, which preserves the 18 outside
substitutions and compensates in the last two cells, reaches the moment
target numerically only with a much larger affine coefficient (about
9.79). It fails the specified amplitude gate; no rate from that
construction is credited here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_daughter_complement_retune_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
```

The receipt `results/early-pole-daughter-complement-retune-2026-11-27.json`
records the selected cells, achieved defects, trial coefficients and initial
cold-mode differences.
