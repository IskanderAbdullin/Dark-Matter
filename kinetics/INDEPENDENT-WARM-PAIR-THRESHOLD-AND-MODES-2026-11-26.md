# Frozen scalar-pair sign threshold and thermal exchange modes

## A lower sufficient threshold than one half

The [frozen pair inequalities](INDEPENDENT-WARM-FROZEN-PAIR-BARRIER-2026-11-26.md)
give a represented-partner entropy gap g and a positive exterior
upper bound P. Their reported inward net lower floor L_net at
tagged occupation 1/2 is (1−exp(−g)) times a pure inward loss
lower bound L. Holding all other modes fixed, all gain
terms scale as 2(1−f) while the forward loss scales as 2f.
Consequently a sufficient occupation threshold for a nonpositive
scalar-pair rate at every represented tag momentum is

    r = exp(−g) + P/L = exp(−g) + (P/L_net)(1−exp(−g))
    f_threshold = r/(1+r).

The worst threshold across the published frozen states is
2.569e-5 on q32 and 2.568e-5 on q64. For the more conservative
explicit time surrogate it is 0.001002 and 0.004711, respectively;
those surrogate thresholds cannot be used to infer the sign
at the much lower actual cold occupation between steps. At the
q32 7.394049 MeV endpoint the frozen threshold is 1.935e-5,
while the first light cell has minimum occupation 4.0977e-5.
At the common-chart q64 7.754805 MeV replay the frozen threshold
is 2.241e-5 and the corresponding minimum is 4.2887e-5.
Entire shells above this conservative threshold contain at least
11.56% and 11.89% of the respective represented light inventories.
This is not a bound on the integrated sign over the rest of the grid,
on the prompt/elastic terms, or on the time-dependent occupation.

## Physical linear modes at the synthetic thermal anchor

Use the [pair number/heat response](INDEPENDENT-WARM-COMMON-PAIR-NUMBER-HEAT-2026-11-26.md)
at T=17 MeV and μ=−30 MeV and evaluate each bath's full
number/heat susceptibility C independently. The susceptibility
has event vector u=(1,(E−μ)/T) and weight proportional to
f(1+f) for the scalar and f(1−f) for light. The pair dynamics
conserve the joint number and heat and relax relative affinities
with matrix (C_scalar⁻¹+C_light⁻¹) J. The symmetric similarity
problem and the full four-moment dynamics reproduce its two positive
relaxation rates and the two conservation zero modes.

Taking the visible bath also at 17 MeV solely for the Hubble
normalization, the slow and fast selected rates are 3385.87 H and
10137.01 H at the finer rule. Suppressing the number–heat cross
coupling would instead give 5118.23 H and 172022.52 H; the latter
is not an actual decoupled heat eigenmode. These are tangent
rates for two internally thermal reservoirs in one scalar-pair
channel, not mode lifetimes for the nonthermal expanding runs.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_pair_frozen_threshold_reference.py \
  results/warm-pair-frozen-barrier-tracks-2026-11-26.json.gz \
  results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
  results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz
python -W error warm_common_pair_thermal_relaxation_reference.py EOS.dat \
  results/warm-common-pair-number-heat-onsager-2026-11-26.json
~~~

[Frozen and surrogate thresholds](results/warm-pair-frozen-threshold-2026-11-26.json);
[thermal susceptibility and eigenmode diagnostics](results/warm-common-pair-thermal-modes-2026-11-26.json).
