# Probe whether the late selected transient is directly sensitive to the outer pad

The [cut-aligned time test](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
uses a 1.05 GeV daughter momentum pad. Freeze the identical
phase at the .68 GeV start and at Hdt=1e-7 and 1e-6,
and at the .80 GeV start. Compare two stationary parent
folds: one disallows outgoing states past the pad, and
one treats incoming exterior as vacuum and allows outgoing
scalar/light states at any kinematically allowed momentum.
The open comparison changes the parent blocking balance
as well as the angular decay destinations.

All four frozen signed cut fluxes preserve their respective
signs. The relative change is below 1.4e-12 in these
probes; this ratio is for the signed excess at the original
.71482 GeV energy cut. On the open branch, the combined
scalar plus light energy newly emitted past the pad is at
most 2.81e-22 in units of retained dark energy times H for
the sampled phases. It then falls to about 6.28e-26 on
the .68 GeV terminal sample. A second independent fold
of convex excess above **both** pad energies agrees with
its direct outgoing interpretation and is smaller than
the full energy that leaves the pad.

This is a direct frozen check, not re-evolution on an enlarged
mesh. It does not exclude emissions at unsampled times,
outgoing particles subsequently re-entering, or hard
population not present in the selected inputs.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_prompt_pad_opening_audit_reference.py EOS.dat \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json \
  results/warm-late-cut-aligned-hard068-test-step0-2026-11-26.json.gz \
  results/warm-late-aligned-hard068-fine-test-step5-2026-11-26.json.gz \
  results/warm-late-aligned-hard068-terminal-6-2026-11-26.json.gz \
  results/warm-late-aligned-hard080-test-step0-2026-11-26.json.gz
~~~

[Three quadrature orders at the four states](results/warm-late-prompt-pad-opening-audit-2026-11-26.json).
