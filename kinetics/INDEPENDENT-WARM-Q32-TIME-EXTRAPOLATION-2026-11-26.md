# First-order extrapolation of the tested q32 moments remains admissible

Use the [coupled four-level time ladder](INDEPENDENT-WARM-Q32-MATCHED-TIME-LADDER-2026-11-26.md)
at identical finite q32 support to check a numerical
diagnostic, (U_R(h)=2U(h/2)-U(h)). Combine each shell's
comoving number and kinetic energy, not the reconstructed
entropy coefficients. Use the same physical endpoint edges
to recover a positive two-moment Bose or Fermi profile in
every case. Keep newly recorded exterior counts only as
algebraic summed counts, not as an extrapolated population
with a specified energy or time of entry.

Between the two finest pairs, the unextrapolated scalar/light
shell-number changes are (1.35	imes10^{-5}) and
(7.79	imes10^{-6}); the corresponding extrapolated
changes are (7.28	imes10^{-7}) and (2.77	imes10^{-7}),
18.5 and 28.1 times smaller. Relative kinetic-energy changes
using the same normalization improve by factors 31–33.
All three extrapolated moment sets admit the required
positive/Pauli-limited reconstruction, with moment recovery
within (3.5	imes10^{-14}) and combined number defect below
(8.5	imes10^{-14}).

This is consistent with cancellation of a leading
first-order numerical term on the tested ladder.
It is not an accepted time integration, an absolute error
bound, or evidence that extrapolation remains realizable
for other initial data, larger steps, or active exterior
states. Do not use the inferred moments to replace an exact
portable checkpoint.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_q32_time_extrapolation_diagnostic_reference.py
~~~

[Diagnostics](results/warm-q32-time-extrapolation-diagnostic-2026-11-26.json).

### Verification with a separately evolved 32-block endpoint

With the fifth actual endpoint, the finest two extrapolated
pairs differ by (2.11	imes10^{-7}) in scalar count
and (8.47	imes10^{-8}) in light count. These are 32.2 and
46.1 times smaller than the corresponding raw step changes.
The reduction is roughly 54–57 for kinetic energy. Moment
admissibility and the combined count check still pass.
This supports leading-error cancellation on the tested
finite grid only; it does not establish a generalized
extrapolation solver.

Add `--include-32` to the command for the
[five-level extrapolation diagnostics](results/warm-q32-time-extrapolation-with-32-2026-11-26.json).
