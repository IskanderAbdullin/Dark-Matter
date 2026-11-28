# Reproduce the entire .710 six-substep clock with a stabilized entropy check

The [previous local refinement](INDEPENDENT-WARM-LATE-NEAR-CUT-REFINEMENT-2026-11-26.md)
described an early root rejected by the original binary64
entropy check even though independently evaluated entropy
increased. Keep that test and rejection unchanged. In a
separate replay of the same initial state, require entropy
to increase by both stabilized long-double summation and
direct integration of the entropy change along the path
between states. The two successive path rules must resolve
the same positive change above the measured numerical
resolution. This changes the entropy calculation used in
the certification, not the moment equations or the choice
of the time steps.

All 36 backward Euler substeps at the prescribed clock
pass this test. Exactly one is classified as negative by the
old check: the fifth root of the first interval, as previously
diagnosed. The greatest recorded relative number or energy
defect is 2.23e−15. The minimum positive stabilized entropy
change is 1.0572e−22, more than 1.15e5 times the declared
long-double subtraction resolution on this system. This
test requires long-double arithmetic with extra precision.

At the last common time Hdt=1e−7, W/rho is 5.67038e−17,
compared with 5.61544e−17 on the three-substep path (+0.978%).
The endpoint signed current is −6.60236e−11 in units of
initial energy times H, compared with −3.35905e−11 before.
The first sign-change bracket on the **native** six-substep
samples is [7.6666667e−8, 8.8333333e−8], the same sampled
bracket as the separate final-interval refinement. The
endpoint signs have independent angular checks. The result
does not bound residual time-step error, establish the exact
zero time or open the momentum pad.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Extract 710-3, including its original summary, as in the count note.
python -W error warm_late_near_cut_entropy_stable_replay_reference.py 710-3 replay-stable
~~~

[All 36 acceptance records](results/warm-late-near-cut-710-stable-six-summary-2026-11-26.json);
[hashed snapshots at the six interval ends](results/warm-late-near-cut-710-stable-six-snapshots-2026-11-26.json.gz).
