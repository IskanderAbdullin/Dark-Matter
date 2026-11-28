# Two complete four-level charge ladders: count alignment does not align the spectra

The independent third native grid is now evolved from its own empty-daughter start with eight charge subcells per archived charged cell. The scalar chart still has 97 bins, compared with 67 on the fine family; D8 has 2496 charged daughter and 2816 parent bins on the third chart. The final restart clock is matched to the previous D1/D2/D4 ladders and the [fine D8 result](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-EIGHT-CHARGE-SUBCELLS-2026-11-28.md). No resampling across the two native histories was used.

| x | Third D8 scalar yield | Third D8 vs third D4 | Third D8 vs fine D8 |
| --- | ---: | ---: | ---: |
| .01 | 3.47247314e−7 | −0.00336% | −0.04841% |
| .02 | 1.01945185e−5 | −0.04274% | +0.04319% |
| .04 | 1.68687960e−4 | −0.48928% | −0.16759% |
| .05 | 3.54942919e−4 | −0.78614% | −0.27162% |

At x=.05 the two D8 counts differ by only 0.272%, less than D4's 0.864%. But at the same fixed pressure mass the D8 third pressure numerator is still **4.371%** below fine, or **4.113% per scalar** after removing the count difference. The pressure discrepancy changes very little across all four levels. The fitted conditional cold spectra give W1 transport of 2.40945e−7 and total variation 0.376807 at D8, versus 2.38533e−7 and 0.375105 at D4. Reconstructed shapes do not get closer simply because two terminal counts do.

Cold conditioning uses q≤1e−5, a shared exact cell edge. The per-bin counts alone give a conservative W1 support enclosure from **0** to 1.54382e−6 at D8. Some support intervals touch without a strict common interior; we do **not** exhibit a common atomic distribution. The shared-edge mass difference gives no meaningful lower total variation at the recorded precision. Thus the reported 37.68% is for the two piecewise-constant numerical histograms, not a guaranteed shape distance for arbitrary within-cell functions. The pressure result similarly conditions on those numerical reconstructions.

The third D8 late fixed-start time drift in scalar number is 3.26825e−5 relative; all four paired-time checkpoint ledgers pass. It does not control spatial resolution. At k=.0015, using the third D8 two-moment optical enclosure and fixed fitted scalar mass, the stored cold-region worst parent threshold reaches 31.29 in numerical energy units (archived cap 32.00), versus 29.95 for fine D8. This is a kinematic warning about the retained parent tail, not a reconstructed rate or a physical bound.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_peer_charge_shape_transport
# Startup: one 2048-step control through .0075.
python -W error early_pole_cell_endpoint_startup_pilot.py third --div 8 --steps 2048 \
  --one --fast-capacity --state-path /tmp/third-d8-start.json
# Continue matched intervals with --restart previous-state, --x-end target,
# --steps 512 then 1024 each, omit --one and keep --div 8 --fast-capacity.
p=results/early-pole-cell-endpoint; d=2026-11-28
python -W error early_pole_endpoint_paired_charge_ladders.py \
  "$p-fine-charge-ladder-d8-$d.json" "$p-third-charge-ladder-d8-$d.json"
python -W error early_pole_endpoint_peer_charge_shape_transport.py \
  --fine "$p-fine"{,-d2,-d4,-d8}"-x05-state-$d.json.gz" \
  --third "$p-third"{,-d2,-d4,-d8}"-x05-state-$d.json.gz"
```
