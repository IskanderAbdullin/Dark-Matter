# Bound the exclusive neutral pion and photon within the weighted count budget

The [weighted direct count](INDEPENDENT-WARM-DIRECT-THREE-PION-DISPERSION-2026-11-26.md)
with charged-pion channels can be extended by one
more distinct vacuum channel, π⁰γ. Its massless
incoming photon has no fixed positive energy floor lab-side.
At fixed invariant mass squared s and total lab energy Q it
instead satisfies E_gamma≥(s−m_pi0²)/(4Q). The photon
Bose factor divided by its MB weight is bounded by
1+4TQ/(s−m_pi0²), using 1/(1−exp(−x))≤1+1/x.
Bound the pion factor using its mass, and integrate
the additional factor Q as an exact Bessel energy
moment. Thus no fictitious photon mass regulates the
soft enhancement. At threshold the vacuum form factor
gives R proportional to (s−m_pi0²)^3/s². Under a
stated uniform form-factor cap, its threshold
suppression absorbs the inverse soft factor.

The [Hoid–Hoferichter–Kubis fit](https://arxiv.org/html/2007.12696)
gives the positive π⁰γ dispersion contribution up to
1.35 GeV. Use a conditional cap of 45.0×10⁻¹¹, the
reported central result plus twice the total quoted
error. A strictly positive subset can use this whole
cap, without asserting the cross section outside
its experimental coverage. The comparison treats
three alternative splits as alternatives, never sums them.

For a separate three-band comparison, use the explicit
assumption |Fhat|≤10 GeV⁻¹ below .32 GeV in the
[prior form-factor convention](INDEPENDENT-WARM-NEUTRAL-PION-PHOTON-BUDGET-2026-11-26.md).
For .32–.60 GeV use a conditional positive moment cap
of 2.6×10⁻¹¹ (twice the fit's sub-.60 GeV result).
For .60–1.35 GeV spend at most the full 45.0×10⁻¹¹.
These transfer windows are disjoint. The caps are
explicit stress choices, not statistical intervals or a
measurement of the threshold form factor.

| Reheating (MeV) | Prior weighted upper | Added neutral-channel upper, three-band option | Updated conditional upper |
|---:|---:|---:|---:|
| 50 | .00076306 | .00001209 | .00077515 |
| 75 | .00428472 | .00005492 | .00433964 |
| 100 | .01514855 | .00014827 | .01529682 |

The 100 MeV addition is less than 1% of the
preceding upper; the full conditional ceiling is
about 25.5% of the .060 half-dark comparison.
There is a separately tabulated 1.35–4 GeV weight
per unit spectral R, but it is not included in the
table. No in-medium spectral function, other
radiative or continuum channels, or dark heavy
regeneration are bounded here.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_neutral_pion_photon_count_reference.py EOS.late.dat \
  "$D/warm-direct-vacuum-temperature-bracket-2026-11-26.json" \
  "$D/warm-direct-three-pion-count-2026-11-26.json"
~~~

[Machine-readable cases, split choices and resolution controls](results/warm-direct-neutral-count-2026-11-26.json).
