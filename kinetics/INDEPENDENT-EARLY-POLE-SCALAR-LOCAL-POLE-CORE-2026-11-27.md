# The signed local vacuum pole remainder depends on how the archive edge is resolved

This is a numerical sensitivity test in the same frozen tagged model as the [guarded annuli](INDEPENDENT-EARLY-POLE-ARCHIVE-SCALAR-NEAR-POLE-2026-11-27.md). It does not add a new term to the live solver. Use the dimensionless invariant offset `x=s/M²−1`, threshold `|x|<R`, and a constant illustrative width `gamma` in x. At fixed x independently fold the squared-pole numerator S(x), its real interference coefficient I(x), and the crossed square U(x). Subtract only the central constant S(0) over the same symmetric interval. The remainder is the integral on x ∈ [0,R] of

```text
 [S(x) + S(−x) − 2 S(0) + x (I(x) − I(−x))] / (x² + gamma²)
 + U(x) + U(−x).
```

This is a signed remainder under one declared convention, not a positive scattering probability. In particular it does not complete the real-particle kinetic matching discussed in the [factorized pole test](INDEPENDENT-EARLY-POLE-SCALAR-POLE-FACTORIZATION-2026-11-27.md). Samples with gamma/R = 0.01 and 0.1 are numerical sensitivity controls, not widths inferred for the medium. A manufactured polynomial gives an independent analytic integral for each term. Invariant quadrature is segmented at the incoming shell edge thresholds on both sides.

At tagged q=4.1e−6, compare the original archive steps with a localized ramp in p³ near one daughter edge (p=1.422222… on fine and p=1.44 on third). This ramp changes only the two adjacent cells and preserves their **combined number and energy** using two small shoulders. It is represented by a finite number of subcells and preserves Pauli-admissible occupations. The table reports the exterior-partner signed remainder at gamma=0 in the reference H normalization.

| Archive and local shape | R=2e−4 | R=2e−3 |
| --- | ---: | ---: |
| Fine, original steps | +1.156e−19 | −5.208e−17 |
| Fine, conserved ramp, eight subcells | −4.196e−18 | −1.025e−16 |
| Fine, conserved ramp, sixteen subcells | −2.958e−17 | −1.315e−16 |
| Third, original steps | +2.332e−19 | +4.141e−18 |
| Third, conserved ramp, eight subcells | +1.816e−19 | +1.011e−17 |
| Third, conserved ramp, sixteen subcells | +1.780e−19 | +9.967e−18 |

For the archive steps at the widest window, increasing the quadrature orders leaves the opposite signs unchanged. On fine, rounding one edge flips the narrow-window sign and the two ramp resolutions disagree substantially in magnitude. Matching these two moments does not identify this tiny remainder or a continuum value. Unresolved medium interference, the actual retarded width and a consistent treatment with an explicit evolving parent are not represented.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
 python -W error early_pole_archive_scalar_pole_core_reference.py "$g"
 for n in 8 16; do
  python -W error early_pole_archive_scalar_pole_core_reference.py "$g" \
   --radii .0002 .002 --width-fractions 0 --radial-orders 12 18 --core-orders 8 \
   --ramp-fraction .03 --subcells "$n"
 done
done
```
Raw: `results/early-pole-archive-scalar-local-pole-core-2026-11-27.json.gz`.
