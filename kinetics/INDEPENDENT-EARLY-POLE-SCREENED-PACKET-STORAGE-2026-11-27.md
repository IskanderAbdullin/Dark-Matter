# Smaller screened packet storage without changing the accepted graph

Unlike and same-flavor elastic geometry now store checked cell indices as 32-bit integers and finish one exchange-momentum batch at a time. Stored weights and all continuously valued coordinates stay 64-bit. The fourth outgoing leg is checked as before when assembling the graph but its unneeded raw arrays are released. The unlike raw momenta needed by the active polarization and remapping diagnostics are still retained.

For both native archives, an independent comparison with the preceding implementation prepared identical packet counts, order, coordinates, cell indices and weights. A SHA-256 fingerprint of every canonical array agrees bit-for-bit after normalizing the storage width of the indices. The comparison includes the unlike raw parent and daughter momenta. An illustrative fine construction with three targeted subdivisions and 1,903,869 packets has 133,714,836 visible array bytes instead of 156,561,264 (14.6% lower). On one Linux cloud run the full constructor peak moved from 729,284 to 708,584 KiB (2.8% lower). An unbatched intermediate with the same final array bytes peaked higher, at 753,276 KiB; final stored bytes alone do not establish peak memory. These are construction fixtures, not an end-to-end speed or maximum-mesh claim.

Reproduce from `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
for g in fine third; do
  python early_pole_packet_storage_check_reference.py "$g" \
    --expected results/early-pole-screened-packet-storage-2026-11-27.json
done
python early_pole_packet_storage_check_reference.py fine --subdivision 3
```

Raw comparisons and construction memory observations: `results/early-pole-screened-packet-storage-2026-11-27.json`.
