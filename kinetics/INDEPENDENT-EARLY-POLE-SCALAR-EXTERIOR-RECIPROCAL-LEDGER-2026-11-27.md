# Exterior decay packets close a reciprocal hard and scalar ledger

For the [matched final five-channel snapshots](INDEPENDENT-EARLY-POLE-JOINED-SCALAR-EXTERIOR-SOURCE-2026-11-27.md), deposit each exterior decay packet into the participating parent and daughter cell's number and energy moments. This is still a single frozen derivative, not a trajectory with the exterior switched on. The elastic screening, hard support and scalar pole are unchanged. The scalar source is recomputed after imposing both triangle inequalities, not imposed from the saved number-only result.

At q=[.1,5] the native-step model gives hard parent energy −6.354e−7, daughter energy +3.664e−7 and scalar energy +2.690e−7 in the reference rate normalization. The scalar receives about 42.3% of the energy lost by the parent in this band. The continuous archive reconstruction gives −6.403e−7, +3.691e−7 and +2.712e−7, respectively. In the lower exterior interval [1e−5,.001] only about .06% of the parent energy loss appears in the scalar. These percentages refer to this decay contribution, not to the total parent energy of the box.

Over the three exterior bands and both reconstructions, the per-species particle ledger closes to 1.1e−22 absolute or better and the reciprocal energy ledger to 1.4e−22, with positive FD/BE entropy production on every band. These are numerical packet identities, not continuum precision. The independent number and energy sums reproduce the earlier exterior source panel within floating-point roundoff. This gives a pathway for a future conservative warm-tail extension; no exterior feedback or model abundance is determined here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for p in steps cell-conservative; do
  python -W error early_pole_scalar_exterior_reciprocal_ledger_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz --projection "$p"
done
```
Raw: `results/early-pole-joined-fine-216-scalar-exterior-ledger-2026-11-27.json`.
