# Eight cold subcells reduce the observed momentum transport without resolving the spectrum

Propagate an independent S8 fine scalar history on the same fixed two-subcell charge chart. It has 249 scalar bins and keeps the original 448 daughter and 488 parent nodes. Every native scalar cell fully below q=1e−5 is split into eight parts before the empty-daughter startup. The startup uses the same 2048/4096 time ladder as S2 and S4, and all later restarts use the same time step counts and retained finer states. This is a comparison of separate evolutions, not a remap of S4.

| x | S8 scalar yield | S8 relative to S4 | TV on common native bins |
| --- | ---: | ---: | ---: |
| .01 | 3.47495103e−7 | −1.098e−5 | .003037 |
| .02 | 1.02181427e−5 | −9.103e−7 | .003048 |
| .04 | 1.73209269e−4 | −2.658e−7 | .003061 |
| .05 | 3.71627215e−4 | −6.230e−7 | .003072 |

The final yield shift remains below the fine-from-fixed-restart temporal difference `3.19e−5`; it does not estimate a continuum yield error. Both number and continuity budgets pass, with final-segment maximum absolute continuity discrepancy 2.78e−11 at retained samples.

On the same [reconstructed histogram support](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-TRANSPORT-AND-NODAL-MOMENTS-2026-11-28.md), the matched S4→S8 result is **6.850%** full-support TV, versus 14.219% for S2→S4. The maximum full-support CDF gap decreases from 4.437 to 1.338 percentage points. The W₁ distance of the conditional cold counts goes from 3.873e−8 to 8.879e−9 (S4→S8 is .245% of the S4 cold mean); the latter mean shifts by −6.630e−9. Thus both measures decrease, but the remaining changes are shorter in q than the histogram TV alone suggests. These are finite-grid contrasts and not an asymptotic bound. The cold number spectrum remains unresolved.

The fixed-mass pressure functional changes by only −5.26e−6 relative, with q≥.01 still carrying 99.8558% of it. This is not evidence that the pressure-contributing warm grid has converged: the S8 experiment leaves that grid unchanged. A [separate warm-only static fixture](INDEPENDENT-EARLY-POLE-ENDPOINT-WARM-NODAL-FIXTURE-2026-11-28.md) tests the nodal representation while keeping the entire histogram fixed.

Replay with `early_pole_endpoint_scalar_subshell_pilot.py fine --charge-div 2 --scalar-div 8`, using 2048 coarse steps to .0075 and 512, 1024, 1024, 1024 at the later fixed epochs. The optional exact capacity fast path was used. The [validated chain](results/early-pole-fine-charge2-scalar8-chain-2026-11-28.json), [x=.05 coarsened comparison](results/early-pole-fine-charge2-scalar8-ladder-x05-2026-11-28.json), [transport receipt](results/early-pole-endpoint-scalar8-shape-transport-2026-11-28.json) and [moment receipt](results/early-pole-endpoint-scalar8-spectral-moment-budget-2026-11-28.json) include the exact input/exposure metadata. Replay states and both temporal arms are archived. Do not use these diagnostics as a CLASS input certificate or infer a physical primordial distribution.
