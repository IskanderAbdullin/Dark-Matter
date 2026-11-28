# Resume accepted expanding blocks without losing the parent or exterior ledger

Long numerical runs need a restart point beyond the daughter
snapshot. Reconstructing the parent from a new equilibrium or
discarding collisionless outgoing cohorts would change a resumed
run. The [comoving driver](INDEPENDENT-WARM-FINITE-PARENT-COMOVING-EXPANSION-2026-11-26.md)
can now opt in to an atomic accepted-block checkpoint. It stores
the full represented parent occupation, the cumulative expansion
work, positive exterior cohorts and their trace, particle-exchange
ledgers, and the accepted daughter phase fingerprints. The
default one-shot numerical path is unchanged.

The restart key binds the EOS bytes, all locked input-guide bytes,
numerical options, both parent-grid fingerprints, and selected
numerical source files. Resume validates the complete accepted
prefix before continuing. A writeout beyond the accepted prefix
is not a checkpoint and is regenerated if present. Concurrent
cooperating writers in the same checkpoint directory are excluded
by a lock. This is a reproducibility and crash-recovery control,
not an authentication mechanism or a guarantee of identical
floating-point output on different software/hardware.

## Restart comparison on identical .68 GeV inputs

On the selected coupled window, stop after one accepted block.
Resume with the same common chart and exact original clock.
Both accepted daughter-step JSON files and the complete summary
agree byte for byte with the uninterrupted backward-Euler run.
The tested phase has 15,450 heavy-parent nodes and retains eight
positive exterior pair cohorts at the final endpoint. Separate
checks reject a mismatched input lineage and an independently
altered accepted output without rewriting the altered file.

The separate SDIRK2 test also reproduces both step files and the
summary byte for byte. In this test an uncommitted second-step
file is placed in the restart directory before resuming. The
first accepted record is validated unchanged and the uncommitted
placeholder is regenerated from the accepted checkpoint. No
parent equilibrium reinitialization or transfer from a different
chart is used.

~~~sh
# Prepare the shared turn12 inputs for the coupled clock note.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
run() { python -W error warm_finite_parent_comoving_expansion_reference.py EOS.dat "$@"; }
for method in be sdirk2; do
 run "uninterrupted-$method" expand-finite-all turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2 --collision-clock "$method"
 run "restart-$method" expand-finite-all turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2 --collision-clock "$method" --checkpoint --stop-after 1
 if [ "$method" = sdirk2 ]; then
   printf '{"uncommitted":true}\n' >"restart-$method/step-2.json"
 fi
 run "restart-$method" expand-finite-all turn12/step-{8,9,10}.json --order 2 --common-subdivisions 2 --collision-clock "$method" --checkpoint --resume
 python -W error warm_comoving_durable_restart_audit_reference.py "uninterrupted-$method" "restart-$method"
done
~~~

[Accepted output hashes and refusal diagnostics](results/warm-durable-comoving-restart-audit-2026-11-26.json).

## Check the passive-parent approximation at saved physical phases

The original integrations use an initial parent distribution
carried collisionlessly in the Friedmann denominator, even when
the actual reacting parent is different. The new checkpoint
contains the actual occupation, so the effect of replacing that
passive parent density can now be evaluated at the *same* daughter,
radiation and exterior phase. Fold the parent difference per
comoving node before summation. Compute the resulting Hubble
shift from the rationalized square-root difference rather than
subtracting two Hubble rates that may round to the same number.

| Saved phase | Instantaneous (H_active − H_passive) / H_passive |
| --- | ---: |
| .68 GeV backward-Euler window endpoint | −3.10e-17 |
| .68 GeV SDIRK2 window endpoint | −3.16e-17 |
| .80 GeV Hdt = 2.50e-7, fourth accepted block of an early-refined test | −2.35e-16 |

This supports using the same prescribed clock for comparisons at
these sampled states. It is not a rerun with the dynamic parent
density in the expansion loop, a global bound over the entire
trajectory, or a bound on additional visible and dark sources.

~~~sh
python -W error warm_checkpoint_parent_gravity_reference.py EOS.dat turn12/step-8.json turn12/step-10.json restart-be
python -W error warm_checkpoint_parent_gravity_reference.py EOS.dat turn12/step-8.json turn12/step-10.json restart-sdirk2
# With guides080-24 generated from the .80 anchor as for the earlier clock refinement:
g=guides080-24
run early080-sample expand-finite-all "$g"/step-{0,1,2,3,4,6,8,10,12,15,18,21,24}.json --order 2 --common-subdivisions 2 --checkpoint --stop-after 4
python -W error warm_checkpoint_parent_gravity_reference.py EOS.dat "$g/step-0.json" "$g/step-24.json" early080-sample
~~~

[Checkpoint gravity comparison](results/warm-durable-parent-gravity-2026-11-26.json).
