# Robustness of the frozen scalar-pair sign to entropy perturbations

Take the [frozen and explicit interpolant bounds](INDEPENDENT-WARM-FROZEN-PAIR-BARRIER-2026-11-26.md)
and keep the momentum support and vacuum incoming exterior fixed.
Allow arbitrary reconstructed entropy perturbations bounded in
sup norm by a for the scalar and b for the light partner. These
are perturbations of log inverse occupation odds. The tagged
occupation stays at 1/2, and Bose entropy remains positive.

The shellwise blocking gap falls by at most 2a+b. The represented
forward-loss lower bound scales by at least the factor exp(−b)
times the ratio of the new to the old blocking margin,
(1−exp(−gap_new))/(1−exp(−gap_old)). For the positive vacuum-partner
channel, the upper scalar product scales by at most exp(2a)
times the square of (1−exp(−z_min))/(1−exp(−z_min+a)), where z_min
is the original minimum scalar entropy in its kinematic band.
No differentiability assumption on the perturbations is needed.

For a=b=3.5, all frozen q32 and partial q64 states retain a
positive gap of at least 0.0701. The worst positive-exterior to
negative-interior ratio is at most 0.00537 on q32 and 2.33e-14 on
q64. On the more conservative time-surrogate boxes, the guarantee
holds with a=b=2 on q32 and a=b=1.75 on q64. Larger tested boxes
that fail the sufficient bound do not demonstrate a sign reversal.

As a numerical scale check, the accepted q32 native/corrected
endpoint differs by at most 1.11e-6 scalar and 7.68e-3 fermion
entropy units on the common grid. The matched q64 common-chart
quadrature change is below 1.52e-3 and 2.55e-4. These are
discrepancies between particular numerical reconstructions,
not error bounds to the true distribution or the live in-between
state. The result still applies only to the scalar-pair term and
does not cover untracked incoming phase space.

~~~sh
python -W error warm_pair_entropy_perturbation_barrier_reference.py \
  results/warm-pair-frozen-barrier-tracks-2026-11-26.json.gz \
  --compare results/warm-multishell-trajectory-native24-q32-2026-11-26.json.gz results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
  --compare results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz results/warm-pair-common-chart-q64-terminal-fine-2026-11-26.json.gz \
  --compare results/warm-pair-common-chart-q32-terminal-2026-11-26.json.gz results/warm-pair-common-chart-q32-fine-terminal-2026-11-26.json.gz
~~~

[Bounds and discrepancy checks](results/warm-pair-entropy-perturbation-barrier-2026-11-26.json).
