# Hard same-flavor packet assembly limits the size of angular work arrays

The finite central-angle hard Møller/Bhabha constructor now subdivides the incoming pair index into batches before forming angular support tables. It preserves packet ordering, cuts, continuous weights and coordinates; only validated cell indices use a narrower integer type. The fourth outgoing leg is still checked before its transient arrays are released. This is independent of the separate screened Q-window packet packing.

On both archived grids with order (1,3,4,1), every canonical graph lane has the same SHA-256 before and after batching. The initial nonthermal collision moment vector, entropy contribution and energy residual also agree exactly. The two fine graph sizes remain 509,512 / 512,008 events, and the two third graph sizes remain 1,084,100 / 1,086,596. In one Linux construction comparison peak resident memory went from 305,604 to 227,640 KiB on fine, and from 525,068 to 349,604 KiB on third. At the higher order (2,4,5,2) the new constructor also completes a fine graph and independent initial collision ledger. Lower construction memory is not a guarantee of total integration memory or of scientific convergence.

Reproduce from `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
for g in fine third; do
  python early_pole_same_hard_graph_fixture_reference.py "$g" \
    --expected results/early-pole-same-hard-packed-geometry-2026-11-27.json
done
python early_pole_same_hard_graph_fixture_reference.py fine --order 2 4 5 2
```

Raw fingerprints, checks, and measurements: `results/early-pole-same-hard-packed-geometry-2026-11-27.json`.
