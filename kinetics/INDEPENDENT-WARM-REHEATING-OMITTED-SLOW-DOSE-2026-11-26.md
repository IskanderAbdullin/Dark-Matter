# Measure the slow dose during the fast return rather than extrapolating a local rate

Evaluate selected slow channels on ten states of the
[actual fast-return trajectory](INDEPENDENT-WARM-REHEATING-FAST-PROMPT-RETURN-2026-11-26.md)
at the 100→60 MeV anchor. This comparison integrates over
only the first 1.5e−8 Hubble times. Each state comes from the
same start with the same prescribed short integrator and
fixed energy, not from extrapolating the final rate.

Fold gross heavy regeneration from both light and scalar
pairs. Separately compute the signed light-antilight ↔
scalar-pair contribution with full final blocking and
stimulation on the **nonthermal** reconstructed shapes.
Retain the analytical opposite energy ledger after prompt
decay. Heavy regeneration is a gross forward-only
contribution in this diagnostic. Inverse heavy channels
were tested on a separate thermal closure and are not
silently added to this nonthermal fold.

| Integrated comparison, divided by visible entropy | Ten sampled clocks |
|---|---:|
| Gross pair events making heavy from scalars | 1.4563e−11 |
| Gross pair events making heavy from light | 3.2566e−11 |
| Scalar count change from the signed light–scalar pair channel alone | −1.0374e−10 |
| Light equivalent with both included terms | +1.3287e−10 |
| Scalar equivalent with both included terms | −3.8608e−11 |
| Scalar energy change in the same projection | −1.0103e−11 GeV |

The two projected count changes are +1.108e−6 of
the starting light count and −6.577e−7 of the starting
scalar count. The estimates are small on this tiny fast
interval despite the large rates per full Hubble time.
Reduce the time grid from ten to six sample points and
the projected scalar change moves by 1.7%; the separate
momentum/angle grid changes it by less than 0.1% on the
ten chosen clocks. This is a sampled separation check,
not a rigorous bound on the unknown unsampled or omitted
channels. It supports testing reduced descriptions, but
does not prove an invariant slow manifold or justify any
prescribed closure over the whole subsequent evolution.
[Invariant-manifold methods](https://arxiv.org/abs/cond-mat/0311017)
motivate that distinction.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error warm_prompt_reheating_omitted_dose_reference.py EOS.late.dat \
 local-spectra/0.100-0.060-171.json local-fast-samples
~~~

[Clock and fold controls](results/warm-prompt-reheating-omitted-dose-2026-11-26.json).
