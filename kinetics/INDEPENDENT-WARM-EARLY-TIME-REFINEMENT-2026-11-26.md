# Time refinement at fixed q16 support

The [early support comparison](INDEPENDENT-WARM-EARLY-SERIAL-SUPPORT-2026-11-26.md)
compared q16 and q32 on four prescribed outer blocks. This test
holds the **q16** support, initial shape, radiation clock, signed
operators, and Jacobian choice fixed. Double the number of outer
blocks three times and match all runs at the original four
expansion endpoints. This tests time stepping, not momentum
discretization, radiation inputs, or unknown prior tails.

| Outer blocks | Inclusive scalar Ωh² at 10.442 MeV | Inclusive stable Ωh² at 10.442 MeV |
| ---: | ---: | ---: |
| 4 | 0.0001622275220 | 84.36601067 |
| 8 | 0.0001624413069 | 84.31629325 |
| 16 | 0.0001625401869 | 84.29329792 |
| 32 | 0.0001625878612 | 84.28221086 |

The final 16→32 changes are **+0.02933%** in scalar count and
**−0.01315%** in stable count, still about 1.14 times the
four-block q16→q32 changes. The log₂ ratios of successive
time-grid differences range from 0.92 to 1.12 at the matched
probes, consistent with a substantial first-order component
in this configuration. These are apparent slopes from actual
discrete changes, not a proven time-error interval or a
Richardson-extrapolated physical abundance. The slowest
component of the split scheme has not been isolated by this
test. The q32 refined results are a separate comparison;
do not subtract q16 at 32 steps from q32 at four steps and
call the gap a momentum-support effect.

## Reproduction

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
S="$D"/50MeV-from-empty-pair-quantum-shape-2026-11-26.json
for n in 4 8 16 32; do
  python -W error warm_coupled_serial_restart_reference.py EOS.dat "$S" "q16-t$n" \
    --stop 1.15 --steps "$n" --local-prompt-jacobian
done
python -W error warm_serial_time_convergence_reference.py \
  "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json "$S" \
  "$D"/warm-early-serial-support-audit-2026-11-26.json \
  q16-t{4,8,16,32}/checkpoint.json
~~~

The archived compressed checkpoints are deterministic gzip copies
and can be passed directly to the final audit in place of plain JSON.
