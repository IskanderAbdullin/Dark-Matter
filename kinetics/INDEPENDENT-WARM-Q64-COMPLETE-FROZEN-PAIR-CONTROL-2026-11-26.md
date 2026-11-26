# Extend the conditional frozen pair sign test to all corrected q64 stages

The [full wider-support trajectory](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md)
allows extension of the prior frozen pure scalar-pair comparison
to every one of its 25 saved states, and over 24 prescribed
time interpolation brackets. Use the existing angular-independent
entropy bound and exterior geometric allowance at tagged occupation
1/2, with the other phases held fixed and exterior partners vacuum.

On the saved states the smallest represented entropy gap is 10.5701,
and the largest permitted exterior positive term is 4.32e-20 of the
interior lower loss. The least-negative uniform upper derivative is
−1.736e-27 GeV in the convention of the frozen pair test. The
bilinear interpolation test also succeeds in all 24 brackets,
but cannot constrain a real intermediate state that does not
follow the prescribed interpolation. Neither result includes
prompt or elastic contributions, an interacting heavy population
or an arbitrary incoming exterior distribution.

Using `ladder/wide` from the linked support reproduction:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_pair_frozen_barrier_track_reference.py ladder/wide/guide-start.json ladder/wide/q64 24 --style guided
python -W error warm_pair_interpolant_barrier_reference.py ladder/wide/guide-start.json ladder/wide/q64 24
~~~

[Frozen states](results/warm-q64-corrected-complete-frozen-pair-2026-11-26.json.gz);
[prescribed surrogate brackets](results/warm-q64-corrected-complete-pair-brackets-2026-11-26.json.gz).
