# Add the exclusive three-pion vacuum channel with a positive data moment

Increase the [restricted weighted count upper](INDEPENDENT-WARM-DIRECT-RESTRICTED-COUNT-UPPER-2026-11-26.md)
by a separate vacuum π⁺π⁻π⁰ input for both light and
heavy pairs. The charged two-pion and three-pion
channels are distinct and are not fitted or counted
twice. For a nonnegative bare e⁺e⁻ cross section,
multiply its experimental dispersion moment by the
largest upper thermal-kernel to lower muon-kernel
ratio over the same covered energy window. This is
valid without assigning the narrow resonance peaks
to points, extrapolating the cross section, or
assuming their internal profile. In Maxwell–Boltzmann
statistics the product of the three incident weights is
exp(−E_total/T); the Bose–Einstein correction is bounded
by independently using the smallest possible energy,
its rest mass, on each of the three pion legs. The
thermal numerator uses the light kernel plus twice
the heavy kernel to account for the weighted ledger.

Use independent [SND Table 7](https://arxiv.org/html/2603.01635#S7.T7)
and [BABAR data](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.104.112003)
alternatives, not a combined fit. Each quotes disjoint
band moments from .62 to 1.10 GeV and from 1.10 GeV
to 1.975 GeV (SND) or 2.00 GeV (BABAR).
Use the central value plus twice the sum of reported
errors as a conditional cap, rather than a confidence
interval. For SND the printed point cross section in
Eq. (14) differs by a factor of three from the usual
unit R ratio. The calculation allows that entire
factor of three in the conservative direction for
SND, rather than assuming that the formula is a typo.

| Reheating (MeV) | Baseline Ω equivalent upper | With BABAR measured bands | With SND measured bands and factor 3 padding |
|---:|---:|---:|---:|
| 50 | .00074844 | .00075336 | .00076306 |
| 75 | .00387454 | .00401276 | .00428472 |
| 100 | .01287463 | .01364154 | .01514855 |

At 100 MeV the two alternatives bound only the
reported windows. In the more padded SND comparison
a uniform unmeasured low-band R value of about 23.2
would be needed to spend the remaining allowance in
this endpoint upper alone if nothing else changes.
This is not a measurement or a claim of such a cross
section. As a separate stress test assuming R≤1
both from the three-pion threshold to .62 GeV and
above the measured bands to 4 GeV, the upper is
.01707968 with SND padding. This stress test is
truncated at 4 GeV and does not bound the omitted
tail or any other hadronic channel. The weighted-count
recycling exclusions, prescribed visible clock and
vacuum assumption are inherited; no point here
bounds in-medium production or establishes a relic.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_three_pion_moment_bound_reference.py EOS.late.dat \
  "$D/warm-direct-vacuum-temperature-bracket-2026-11-26.json" \
  "$D/warm-direct-restricted-count-upper-2026-11-26.json"
~~~

[Both grid controls, experimental alternatives and all bands](results/warm-direct-three-pion-count-2026-11-26.json).
