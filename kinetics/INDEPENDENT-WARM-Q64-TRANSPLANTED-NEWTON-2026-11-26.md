# An archived narrower shape can seed the same wider nonlinear root

The wider q64 warm continuation develops tiny explicit outer
populations. Its second block, matched to the actual older q32
10.339 to 10.288 MeV interval, needed an alternate nonlinear
initialization before accepting. Test whether the nearby
**already computed** q32 endpoint can serve as a numerical
initial guess without replacing the wider physical equations.

Prepare the exact q64 midpoint: half redshift, first signed
pair kick, same radiation and archived exterior gravity,
selected elastic panels and prompt Jacobian. Change only the
initial entropy coordinates of the shared shells to their
archived q32 end values. Keep all q64 outside-shell guesses,
physical moment targets, quadratures and conservation checks.
After the implicit solve, replay the second pair kick and
the final redshift. Compare the retained moments against the
actual accepted q64 block, which was independently computed
using the original solver.

The borrowed, row/column equilibrated Newton iteration
reaches (4.2	imes10^{-12}) normalized residual in four
updates. Compared with the original accepted wider block,
the retained shell-number (L^1) discrepancies are
(1.05	imes10^{-14}) and (8.65	imes10^{-15}) of their
respective species counts, and total retained energy differs
by (1.02	imes10^{-14}) relative. The wide result is not
a substitution of the narrow distribution.

The initial Jacobian probes alone do not explain runtime:
normalizing the predicted residual by shell numbers exposes
a very broad q64 range, but the worst initial residual for
the direct input is in an existing shared shell. This test
does not identify a general performance mechanism or a
verified resource failure. It validates one reusable
initialization on one measured block, with no change to
the production integrator. No elapsed-time speedup or
same-root guarantee on other histories is asserted.

~~~sh
zcat results/warm-extended-coupled-q64-1-blocks-2026-11-26.json.gz > q64-prev.json
zcat results/warm-extended-coupled-q64-2-blocks-2026-11-26.json.gz > q64-end.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_q64_transplanted_newton_warm_start_reference.py \
 EOS.dat q64-prev.json q64-end.json
~~~

[Diagnostics](results/warm-q64-transplanted-newton-warm-start-2026-11-26.json).
