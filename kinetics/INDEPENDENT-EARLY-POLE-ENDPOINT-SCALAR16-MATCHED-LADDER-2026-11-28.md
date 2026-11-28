# The S16 cold scalar history contracts transport again, but not the warm-grid uncertainty

Extend the [matched scalar-only histories](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR8-MATCHED-LADDER-2026-11-28.md) to 16 subdivisions of every complete cold native shell with upper edge ≤1e−5. The charged chart stays at two subcells; the warm grid is unchanged. The new history uses 457 scalar, 448 daughter, and 488 parent nodes. It starts independently from empty daughters and uses the same 4096-step retained startup and the same later fine-clock counts as the earlier scalar sequence.

| x | S16 scalar yield | S16 relative to S8 | TV on common native bins |
| --- | ---: | ---: | ---: |
| .01 | 3.47494193e−7 | −2.617e−6 | .0007561 |
| .02 | 1.02181404e−5 | −2.298e−7 | .0007589 |
| .04 | 1.73209258e−4 | −6.251e−8 | .0007622 |
| .05 | 3.71627201e−4 | −3.917e−8 | .0007650 |

The final paired frozen-restart time difference is 3.17e−5 relative for total yield, larger than the S8→S16 count change. Both arms and all sampled continuity and pairing tests pass. The combined output cannot resolve a continuum yield correction of the observed size.

On the full reconstructed histogram, TV declines from 14.219% (S2→S4) to 6.850% (S4→S8) to **3.398%** (S8→S16). The maximum cumulative full-support gap in the last step is 0.3454 percentage points. The conditional cold W₁ distance similarly drops from 3.873e−8 to 8.879e−9 to **2.220e−9**, and the last cold mean decreases by 1.677e−9. This is a clear contraction of these particular finite-grid contrasts; it does not establish an asymptotic error bound or prove the actual transfer spectrum.

A frozen forward kernel would multiply each q mode by 1+f. On the exact union of the S8/S16 cold shell interfaces, the absolute relative change of this multiplier exceeds 5% on modes carrying 58.07% of the old cold count, and 10% on modes carrying 23.54%. At the five predeclared cold q probes the changes range from −17.40% to +14.82%. These are count-weighted measures of where a local input changes, not estimates of the integrated or net collision rate. The [multiplier receipt](results/early-pole-endpoint-scalar-bose-multiplier-contrasts-2026-11-28.json) is reproduced with `early_pole_endpoint_scalar_bose_multiplier_contrasts.py` on S2/S4/S8/S16 final states.

The fixed-mass pressure per particle changes by −1.45e−6 relative in the final step while q≥.01 still supplies about 99.856% of that estimator. This is not a test of the warm grid. The [independent warm-only history](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR4-WARM2-EVOLUTION-2026-11-28.md) instead changes reconstructed pressure by percent levels. Read these as separate sensitivity axes, not as one total error estimate.

Replay with `early_pole_endpoint_scalar_subshell_pilot.py fine --charge-div 2 --scalar-div 16 --fast-capacity`, using 2048 coarse steps at startup and 512, 1024, 1024, 1024 in the subsequent joins. The [chain](results/early-pole-fine-charge2-scalar16-chain-2026-11-28.json), [final comparison](results/early-pole-fine-charge2-scalar16-ladder-x05-2026-11-28.json), [transport receipt](results/early-pole-endpoint-scalar16-shape-transport-2026-11-28.json) and [moment budget](results/early-pole-endpoint-scalar16-spectral-moment-budget-2026-11-28.json) preserve complete exact replay and diagnostics. These numerical histograms are not certified CLASS inputs.
