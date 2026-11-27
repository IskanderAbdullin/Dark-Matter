# Exterior scalar gain under live hard channels and distributed shell controls

Continue the [reciprocal exterior experiment](INDEPENDENT-EARLY-POLE-LIVE-SCALAR-EXTERIOR-DECAY-2026-11-27.md) from the same accepted fine hard snapshots. The joined control retains the reciprocal decay packets, hard photon pair, Compton and flavor-conversion terms, unlike screened scattering with a refreshed active response, and held same-flavor polarization. Elastic midpoint solves are surrounded by the matched slow half steps. No expansion is included. This is a conditional short continuation at exposure 2e−6.

| Control, q=.1–5 warm shell selection | Integrated warm gain |
| --- | ---: |
| Native markers, decay only | 5.496615e−13 |
| Native markers, joined n6 | 5.496866e−13 |
| Native markers, joined n12 | 5.496872e−13 |
| Continuous hard reconstruction, joined n6 | 5.522412e−13 |
| Native hard state, distributed warm shell order 4, decay only | 5.228875e−13 |
| Native hard state, distributed warm shell order 8, decay only | 5.228499e−13 |

The joined hard terms increase the native warm gain by about .0046% in this window, and the native n6→n12 comparison shifts the gain by roughly one part per million. The whole high exterior node selection q>1e−5 gives 8.781992e−13 in the native joined n6 test. The independent total hard count and energy ledgers remain within 3e−12 and entropy within 3e−11 in every joined arm. The much smaller difference between the separately advanced cold populations is not resolved by this test.

In a distinct decay-only control, the two arms both use the [energy-matched positive shell reconstruction](INDEPENDENT-EARLY-POLE-WARM-MARKER-SHAPE-2026-11-27.md): either retain those resolved nodes as spectators or bring them live. The shared initial scalar energy agrees with the marker model within 2e−21 while allowing the common gap to change by about +4.46e−14 at fixed coupling. The resulting warm gain is about 4.9% below the single-marker control and changes by less than .008% between warm quadrature orders 4 and 8. Thus the warm source is more sensitive to this initial within-shell representation than to the particular short additional hard evolution. These comparisons do not supply a continuum warm distribution or a combined distributed-and-joined result.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --joined --panels 6
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --joined --panels 12
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --joined --panels 6 --projection cell-conservative
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --joined --panels 6 --outer 1.00001e-5 5
for n in 4 8; do
 python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --warm-order "$n"
done
```
Raw: `results/early-pole-live-scalar-exterior-joined-distributed-2026-11-27.json.gz`.
