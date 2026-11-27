# Exact geometric budget for an unresolved daughter rest strip

For fixed positive daughter masses m,b and prescribed parent excess,
consider a scalar momentum strip 0≤p≤p_c. Keep every allowed fermion
energy; no distribution or occupation factor is folded into this
geometric test. Define D=√(K²−m²b²), with K=(M²−m²−b²)/2 and the
near-threshold factors evaluated from the separately supplied excess.
At fixed scalar momentum the allowed fermion-energy width is 2Dp/m².
With u=asinh(p_c/m) the exact two-energy area is

    A = (D/2)[sinh(2u)−2u].

It is bounded between L/√[1+(p_c/m)²] and L, where
L=(2D/3)(p_c/m)³. Thus the omitted unweighted area at the rest edge
has an explicit cubic bound even when its absolute energy edge
rounds to m. The geometrically weighted scalar trace has an
independent closed form 2D[p_c−m arctan(p_c/m)]. The mean companion
trace is trapped between its endpoint weights; as p_c→0 it tends
to m b²/K.

Compare these formulas to the independent
[kinetic-offset quadrature](INDEPENDENT-TWO-MASS-KINETIC-OFFSET-TRACE-2026-11-27.md)
for two asymmetric sets of masses and five cutoffs from p_c/m=.01
to 10⁻¹⁰. The combined area and scalar trace discrepancy is below
6.32e−11 relative in these tests. At p_c/m=10⁻¹⁰ the scalar
absolute edge has collapsed, but the two unweighted areas are
2.70547e−33 and 1.03575e−32 in the respective fixture units.

This is an exact geometric bound, not a bound on a Bose-enhanced
source by itself. To bound a missing rate, a separate bound on
the absolute collision weight on the strip is required. No bound
on an infrared occupation or on a cosmological relic is assumed.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_rest_strip_exact_area_reference.py
```
[Full receipt](results/two-mass-rest-strip-exact-area-2026-11-27.json).
