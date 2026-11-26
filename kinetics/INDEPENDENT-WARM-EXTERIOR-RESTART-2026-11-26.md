# Carry positive exterior gravity over a warm checkpoint

The earlier [coupled warm benchmark](INDEPENDENT-WARM-COUPLED-B35-2026-11-26.md)
keeps unresolved exterior cohorts during a single stage and audits
their counts across archived stages, but its historical restart files
do not contain those cohorts. A summary count and energy alone do
not determine the later redshifted pressure. This check makes
the *available* positive collisionless cohorts portable without
inventing a distribution for any unarchived tail.

At restart scale A, change each cohort's injection scale to a_i/A,
its comoving count to N_i/A³, and cumulative injected comoving
energy and trace to the old values divided by A⁴. Its rest mass,
mean energy at injection, high/low tag and source channel do not
change. The lower energy uses one mean per cohort. The previously
implemented positive-distribution upper bound is carried and checked,
not used as a second gravitational component. The new warm stage
includes the imported lower energy in the Hubble rate. Its energy
ledger subtracts the entry exterior energy and advances only the
increment in exterior trace. The retained-plus-new-tail count ratio
continues to exclude prior counts, so a multistage count audit can
carry them once, not add them again at each restart.

## Verified handoffs

- An actual *source-only, two-step* portal calculation from 50 to
  45 MeV exports 30 positive cohorts. Their species-resolved counts,
  one-mean energy and upper bound reconstruct the independent source
  ledger. Exporting changes neither the producer report nor the
  original fields of its shape.
- The one-block signed-pair continuation of that shape imports all 30
  cohorts and ends with 36. The normalized stage energy defect is
  about 3.7×10⁻¹⁶ and the total number audit (including prior and new
  tails once) is at the roundoff level.
- A separate 1 MeV check inserts a labeled **synthetic** positive prior
  cohort, then compares two uninterrupted blocks to the same blocks
  separated by serialization and restart. It checks shell number
  and energy moments, the clock, exterior total energy, its upper
  bound, and trace. The runs agree at roundoff, not a
  time-step or phase-space convergence test.

These are controlled numerical handoffs, **not** a rerun of the
full signed 50→12 MeV production with the larger momentum support.
Old tails remain collisionless, do not act as inverse, stimulating,
or re-entering populations, and cannot be reconstructed for older
snapshots that lack cohort data. The default CLI and old shape schema
remain usable unchanged. For an old shape with no cohort data the
new option starts an empty tracking record; it does not silently
replace missing history by inferred particles.

## Reproduction

Use the original EOS described in the
[benchmark setup](INDEPENDENT-BENCHMARK-SETUP.md) and run:

~~~sh
export OPENBLAS_NUM_THREADS=1
common=(python -W error joint_resonant_frw_split_reference.py EOS.dat
  --physical-amplitude --visible-source --visible-backreaction
  --trace-exterior --source-only --empty-pulse 1e-8 --empty-stop .90
  --support-low .005 --empty-cells 18 --matched-vacuum-order --steps 2)
"${common[@]}" --save-snapshots producer-shape.json --save-exterior-cohorts >producer.json
"${common[@]}" --save-snapshots legacy-producer-shape.json >legacy-producer.json
python -W error coupled_warm_pair_continuation_reference.py EOS.dat producer-shape.json \
  --stop 1.0001 --steps 1 --save-snapshots real-bridge-shape.json --track-exterior >real-bridge.json
python -W error warm_exterior_restart_audit_reference.py EOS.dat \
  results/warm-coupled-1MeV-12-and-20-shape-2026-11-26.json \
  producer.json producer-shape.json legacy-producer.json legacy-producer-shape.json \
  real-bridge.json real-bridge-shape.json
~~~
