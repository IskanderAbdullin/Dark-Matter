# Finite-panel cold coefficient bounds conditional on known moments

Use the two verified, complete S2 or S4 panel streams, the actual evolved scalar state, and the original held transverse-pole closure for that chart. Divide the recorded forward coefficient by its fixed within-cell (1+f) and inverse by (f) on nonempty cells to recover the corresponding bare integrated panel weights. Then allow each panel its own nonnegative scalar occupation, while **preserving each original scalar cell's count in the retained band**. In a second calculation additionally preserve each retained cell's mean (E_s-m_s). Original out-of-band values are unchanged, and charges, masses and optical coefficients do not evolve in either test.

A count-only extremizer can place the count in one panel of its cell. With one additional mean constraint, an extremizer uses one or two panels; enumerate the exact two-point weights. To bound the conditional number-neutral occupation, optimize (B-n D) jointly on the **same** feasible occupations and solve each scalar root. Simply dividing independent inverse and slope extrema would give a looser and not generally attainable interval. The recorded baseline lies inside every hull. Both independently computed orders give the same finite optimization endpoints within 3.89e−14 in relative fraction.

| Chart | Gross forward shift, counts only | Gross forward shift, counts and mean excess | Conditional neutral shift, counts and mean excess |
| --- | --- | --- | --- |
| S2 | −2.6643% to +1.2031% | −0.003894% to +0.002506% | −0.001557% to +0.001571% |
| S4 | −0.4976% to +0.5744% | −0.0003481% to +0.0000354% | −0.0000291% to +0.0000437% |

The strong reduction applies when **both moments are prescribed** and only the nonnegative piecewise-constant values on these fixed finite panels can vary. It says nothing about the numerical error in the moments themselves, moving panel boundaries, the unresolved continuous within-panel shape, or a physical variation of the coupled system. In particular it does not reconcile the previously observed local S2→S4 differences between actual trajectories. The narrow conditional root intervals are not predicted photon populations.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_panel_reweight_envelope
d=2026-11-28
for s in 2 4; do
 python -W error early_pole_endpoint_cold_panel_reweight_envelope.py \
   results/fine-charge2-scalar$s-endpoint-continuum-cold-{1,2}.jsonl.gz \
   results/early-pole-fine-charge2-scalar$s-x05-state-$d.json.gz \
   results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json --scalar-div "$s"
done
```
