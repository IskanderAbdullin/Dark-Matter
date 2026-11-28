# An explicit spectator adjustment can restore all three global moments

The [within-shell shape proof](INDEPENDENT-EARLY-POLE-WARM-MARKER-SHAPE-2026-11-27.md) shows that a positive spread around an original marker cannot preserve that individual shell's number, energy and inverse-energy moment simultaneously. Test a different question: after the same warm energy-matched reconstruction, is there a positive change among the **other** archived spectators that restores all three *global* initial moments? Allow 30 unpromoted spectators with q>1e−5, and keep the cold population and the hard state unchanged. This is a manufactured control, not an archive inference.

Minimize the squared relative count changes of those spectators subject to three linear constraints at the original pole: zero added number, no net added free energy, and compensation of the warm inverse-energy surplus. The constraints are row-scaled for numerical conditioning. For both hard reconstructions the nonzero moves have a combined absolute count of about 5.469e−12 and a maximum per-shell fractional change of .0356%; all populations stay nonnegative. The three moment defects are below 1e−24 absolute. Keeping the original coupling then reproduces the initial scalar pole to 1.1e−19 and gives the same evaluated cold source in both controls.

The warm source is still 2.61446e−7 (native hard steps) or 2.62668e−7 (continuous hard reconstruction), about 4.9% below their respective marker sources. Thus the within-shell warm rate difference persists even when the original global initial number, energy and gap are restored through the separate spectator move. That move is not a unique physical prescription and should not be silently applied in an evolved calculation.

From `kinetics/`:
```sh
for p in steps cell-conservative; do
 OPENBLAS_NUM_THREADS=1 python -W error early_pole_scalar_warm_marker_shape_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz --projection "$p" --order 12 --compensate
done
```
Raw: `results/early-pole-warm-global-moment-control-2026-11-27.json`.
