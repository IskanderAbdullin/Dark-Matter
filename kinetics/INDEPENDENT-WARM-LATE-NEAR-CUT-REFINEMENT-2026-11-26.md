# Refine the late .710 GeV reversal from one unchanged saved phase

The [three-substep near-cut trajectory](INDEPENDENT-WARM-LATE-NEAR-CUT-SHAPES-2026-11-26.md)
uses negative endpoint current to bracket a turn only between
Hdt=3e−8 and 1e−7. Keep the accepted phase at 3e−8 byte-for-byte
and refine just that last interval to six equal backward Euler
substeps. The three-substep replay reproduces the originally
saved final snapshot hash. The six-substep replay has positive
current +1.1026e−11 at 7.6667e−8 and negative −3.4591e−11 at
8.8333e−8. Its terminal current is −6.2173e−11, compared with
−3.3591e−11 on the original three substeps. The final W changes
by +0.859% of the original W; the sign of the weak terminal
current is stable in this **local** refinement.

This is a new implicit solve from the same saved state, not
an interpolation or a full-clock convergence check. All six
roots pass the strict entropy guard, conservation, independent
angular projection and the discrete W balance. Neither the
new sampled turning bracket nor the 0.859% change measures
the error from the earlier coarse intervals.

An attempt to use six substeps on every interval was rejected
by the unchanged entropy guard before the first saved positive
time. The final-interval result does not replace that failed
whole-clock test. See the associated early-entropy diagnostic
for the recorded size of the rejected numerical change.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Extract the 710-3 archived folder as in the near-cut count note.
for n in 3 6; do
  python -W error warm_late_near_cut_terminal_refinement_reference.py 710-3 replay-$n.json --substeps "$n"
done
~~~

## Diagnose the rejected early root without accepting it

The fifth of the six proposed roots in the first 1e−11
interval shows a negative change −1.024e−23 in the original
binary64 sum of approximately 5.282e−10. Evaluate the
derivative of each FD/BE entropy along the straight path
in entropy variables, which avoids subtracting near-equal
total values. The 8, 16, and 32 node path rules give the
same positive change +1.05844e−22. An independent stabilized
long-double absolute-entropy difference and the concavity
lower bound agree on the positive sign at the represented
quadrature. The terminal collision entropy times the step
is +1.05784e−22. The binary64 guard had rejected a concavity-
consistent increase. This diagnostic leaves the original
rejection in place; it does not by itself validate the rest
of the six-step-per-interval clock or other platforms.

~~~sh
python -W error warm_late_near_cut_early_entropy_audit_reference.py 710-3
~~~

[Three-substep replay](results/warm-late-near-cut-710-terminal-3-2026-11-26.json);
[six-substep replay](results/warm-late-near-cut-710-terminal-6-2026-11-26.json);
[rejected early-root audit](results/warm-late-near-cut-early-entropy-audit-2026-11-26.json).
