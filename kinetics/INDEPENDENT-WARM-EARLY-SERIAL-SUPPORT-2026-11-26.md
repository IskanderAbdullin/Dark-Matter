# Compare actual q16 and q32 support on the same short warm grid

Start from the two [signed 50-to-12 MeV production states](INDEPENDENT-WARM-SIGNED-Q32-PRODUCTION-2026-11-26.md).
The extended shape retains the same inner edges and five
additional outer shells. Their archived prior exterior tails are
carried **only as counts in the audit**, not supplied to the
subsequent gravitational or collision calculation. The archived
initial shapes lack per-cohort histories. Newly emitted positive
cohorts during this test do contribute to the warm clock, and
are carried through both restart boundaries.

The two runs use the same prescribed log-expansion grid: four
blocks to relative expansion 1.15 (visible T ≈10.442 MeV).
They use the same signed pair kicks, implicit prompt/elastic
operator, parent-local Jacobian and radiation clock. Each run
was stopped at a completed checkpoint and resumed. The result
measures sensitivity to the extended source history and support
together. It does not isolate either effect, test a time-grid
refinement or complete the separate 12→3.43 MeV high-support run.

The audit independently reconstructs the retained number and
newly exterior number from accepted block count changes. At
the endpoint it rechecks both against the stored occupations
and explicitly serialized positive cohorts. Earlier, unresolved
source counts are only added once, with the original mass and
visible-entropy normalization. This does not make those old
tails available for inverse reactions or stimulation.

## Matched endpoint result

| At visible T=10.442 MeV | q16 | q32 | q32 relative to q16 |
| --- | ---: | ---: | ---: |
| Inclusive 86 eV scalar fixed-count Ωh² | 0.000162228 | 0.000162269 | +0.02567% |
| Inclusive stable fixed-count Ωh² | 84.36601 | 84.35633 | −0.01148% |

The stable proxy is still about 703 times the 0.120 reference
at this early point. It is **not a late-time prediction**.
The signs of the support differences persist at all four matched
probes. The last q32 block rejected one initial solver attempt
and accepted the retry before writing the endpoint; the log
does not count the rejected attempt as a completed block.
No other block recorded a retry. Both solver runs stopped by the
same expansion factor and retained the expected absolute phase-space
scales. No prior unresolved shape is reconstructed.

## Reproduction

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
S16="$D"/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
S32="$D"/50MeV-from-empty-pair-quantum-q32-shape-2026-11-26.json
for x in 16 32; do
  if [ "$x" = 16 ]; then s="$S16"; else s="$S32"; fi
  python -W error warm_coupled_serial_restart_reference.py EOS.dat "$s" "q$x" \
    --stop 1.15 --steps 4 --local-prompt-jacobian --max-blocks 2
  python -W error warm_coupled_serial_restart_reference.py EOS.dat "$s" "q$x" \
    --stop 1.15 --steps 4 --local-prompt-jacobian
done
python -W error warm_early_serial_support_audit_reference.py \
  "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json \
  "$S16" q16/checkpoint.json "$S32" q32/checkpoint.json
~~~
