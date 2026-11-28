# Resolving warm marker shells changes the initial inverse-energy moment

The [live warm marker control](INDEPENDENT-EARLY-POLE-LIVE-SCALAR-EXTERIOR-DECAY-2026-11-27.md) preserves occupied native scalar marker nodes exactly. Check what changes if the same eleven q≥.1 shells are instead resolved with positive within-cell quadrature. Keep the accepted final hard snapshot, cold counts, other spectators and Hartree coupling fixed. Compare a uniform occupation in the native q³ measure with a positive exponential reweighting in scalar energy, constrained to retain each marker shell's free-energy moment at the original pole.

For an original marker at energy E0, the identity
`1/E − 1/E0 + (E−E0)/E0² = (E−E0)²/(E E0²)`
shows why keeping both count and mean E exactly cannot also preserve the inverse-energy moment of a point marker when the shell is spread positively. The positive variance expression is independently checked in each run.

| Native fine reconstruction | Warm number rate | Change in scalar pole | Change in cold source |
| --- | ---: | ---: | ---: |
| Archived single nodes | 2.74831e−7 | 0 | 0 |
| Uniform in q³ | 2.61446e−7 | −1.0983e−13 | +2.7521e−11 |
| Energy-matched distributed shell | 2.61446e−7 | +4.4611e−14 | −1.2568e−11 |

Both continuous warm-shell choices lower the frozen warm source about 4.9%; the count-conserving continuous hard reconstruction behaves similarly. Uniform rounding raises the selected warm free energy by 6.9931e−9 (about 3.9% of the marker warm energy); the energy-matched control instead has a first energy difference below 2e−21, including the equilibrated Hartree energy within roundoff. Its warm inverse-energy moment rises by 3.0510e−9 (about 1.1% of the original warm inverse-energy moment). It is incorrect to hold the original scalar gap exactly while also describing this positive within-shell reweighting as an energy-and-number matched reconstruction at unchanged coupling.

Doubling the warm quadrature order 6→12 moves the warm number source by under 0.01% relative. The observed cold-rate shifts are conditional frozen readouts, not a measure of the true late-time feedback or uncertainty of the physical production rate.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
 for n in 6 12; do
  python -W error early_pole_scalar_warm_marker_shape_reference.py "$b" --projection "$p" --order "$n"
 done
done
```
Raw: `results/early-pole-warm-marker-shape-2026-11-27.json.gz`.
