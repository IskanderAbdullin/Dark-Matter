# Add a hypothetical bare scalar mass to the frozen stress, not to the trajectory

Use the same terminal scalar occupations as the
[conformal stress audit](INDEPENDENT-EARLY-POLE-DAUGHTER-RETAINED-KINETIC-STRESS-2026-11-27.md),
including its prescribed exterior spectators, while varying a
hypothetical bare mass (b) at fixed scalar coupling. The
effective gap obeys (m^2=b^2+(lambda/2)sum n/E), with
counterterm (V=(m^2-b^2)^2/(2lambda)).
Then (ho=sum nE-V) and (P=sum np^2/(3E)+V)
give (ho-3P=b^2sum n/E).
An independent dilation check keeps (b) fixed, re-solves
the scalar pole and verifies (P) at each chosen threshold.
At fixed occupations an increasing bare mass raises the
effective gap and energy while lowering both kinetic and
interaction pressure.

The required bare mass in units of the initial massless
effective gap is approximately .667 at (P/ho=.2),
1.132 at .1, 3.017 at .01, and 7.957 at .001. These values
agree across both hard resolutions and the new angular
mesh in the frozen tests. They are *not* an assumed physical
bare mass or the time of a cosmological transition.

Check a different question independently: whether a scalar
with the altered effective pole can still fit the **original
native hard packet geometry**, holding that earlier geometry
fixed. For each selected cold momentum, maximize the parent
minus daughter energy over the allowed relative angle;
the minimum daughter momentum given the parent fixes
the maximizing angle. A convexity check on each smooth
affine-energy branch reduces the calculation to endpoints.
Independent packet searches just below and above the ceiling
verify the support classification.

| Support retained in the frozen hard model | Last effective pole (m/T) | (P/ho) at the ceiling, frozen terminal scalar state |
|---|---:|---:|
| Same selected daughter interval | .00089079–.00089112 | .30640–.30684 |
| Any represented daughter interval | .00277066–.00277099 | .009098–.009101 |

Thus the particular local production interval in the earlier
benchmark loses its accessible cold modes well before the
frozen scalar stress looks matter-like. The last **represented**
hard support can persist to (P/ho) about .0091.
These are conditional comparisons of terminal scalar stress
and earlier fixed native geometry, not a computation of
production with a new bare mass. Unrepresented momenta,
kinetic broadening, other reactions and actual evolution
are not excluded.

Use `early_pole_daughter_bare_scalar_stress_control.py` with
the same three paired controls as the conformal stress
audit. An independent zero-momentum fixture checks the
fixed-bare-mass trace identity.
