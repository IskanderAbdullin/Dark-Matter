# Weight the Bose gain bound by the archived parent occupation

Keep exactly the isolated forward-only assumptions of the
[bootstrapped peak bound](INDEPENDENT-WARM-PAIR-BOOTSTRAPPED-BOSE-CEILING-2026-11-26.md).
The parents may redshift and disappear selectively, but receive no
inverse gains or redistributed occupation. Thus the per-mode parent
occupation at physical momentum p and expansion b is at most
the initial reconstructed occupation at momentum bp. This comparison
includes produced scalars outside the archive, but starts without
unrecorded inherited scalar occupation.

Let `I = integral p f_chi(p) dp` at the anchor. In the integral
for a fixed outgoing scalar, average the first parent's inverse
distance factor over its isotropic directions. Its average is
`1/max(p,k) <= 1/p`. The first parent's radial gain factor is
then at most `4 pi I/(m b²)`. The second parent's angle fixes
the energy delta, and its radial factor is at most `I/(m b²)`.
Replace `P^4 S^2` by `4 I^2` in the earlier gain coefficient;
leave the amplitude cap, outgoing-energy floor, spin/label
overcount and photon-clock assumptions unchanged.

In the archived two-moment reconstruction within a shell,
(p\,dp=(m+K)\,dK) and the Fermi occupation is monotone in the
shell kinetic coordinate (K). Sum conservative endpoint maxima
and minima on increasingly fine partitions. The independent
positive radial quadrature lies within all three nested bounds.

| Anchor | Uniform peak dose | Parent-weighted upper dose | Revised cap for all scalar modes |
| --- | ---: | ---: | ---: |
| 1.206 MeV | 0.001399 | about 1.3×10⁻⁵ | 0.000719 |
| 1.000 MeV | 0.000421 | about 4.7×10⁻⁶ | 0.000716 |
| 0.414 MeV, base | 8.73×10⁻⁶ | about 8.1×10⁻⁸ | 0.000706 |
| 0.414 MeV, padded | 1.83 (no closure) | about 1.3×10⁻⁵ | 0.000719 |

The broad padded support by itself made the previous uniform bound
inconclusive. Accounting for the tiny occupation at large incident
momentum closes the same gain comparison. This does not introduce
a cutoff in the outgoing scalar modes or guess any older exterior
occupation. It remains restricted to ballistic loss-only parents,
not the full signed warm evolution with redistribution and sources.

~~~sh
D=results
OPENBLAS_NUM_THREADS=1 python -W error warm_pair_parent_weighted_bose_ceiling_reference.py \
  "$D"/warm-coupled-pair-b10-28-shape-2026-11-26.json "$D"/warm-coupled-1MeV-12-and-20-shape-2026-11-26.json \
  "$D"/warm-late-paironly-{base,pad}-shape-2026-11-26.json
~~~
