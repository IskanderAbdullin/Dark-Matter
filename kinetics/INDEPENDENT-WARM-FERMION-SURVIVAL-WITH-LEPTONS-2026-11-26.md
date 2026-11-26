# Tighten the warm fermion floor and allow loss to vacuum charged leptons

The [first fermion lower comparison](INDEPENDENT-WARM-ISOLATED-FERMION-SURVIVAL-2026-11-26.md)
used a loose angle-independent vacuum annihilation ceiling and
switched off the visible sector. Keep the same post-10.442 MeV
matched warm reservoirs and the same **imposed** outgoing scalar
occupation cap `F_B`. The fermion ledger counts light and heavy
particles of both charges. Redistribution, transitions between
light and heavy fermions, known tail reentry, and inverse gains
are allowed as before. Count all surviving fermions at least at
the light mass.

## Angular rate bound, distinct from the number-growth coefficient

Enclose the exact crossed vacuum scalar-pair amplitude as a
positive polynomial over its propagator denominators. Use the
[same angle/invariant method](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md)
as in the energy-weighted inverse proof, now directly for
forward `σ v_CM` with either species initially. On each finite
cell use the monotone endpoint factors and the exact
maximum of z(1−z). Above s=2 GeV² the geometric Cauchy
denominator is increasing, giving an analytic decreasing
tail bound. Quotient from CM to lab never increases the rate.
For the light and heavy initial states respectively, the
resulting all-momentum upper coefficients are
2.70205×10⁻¹² and 1.33211×10⁻¹² GeV⁻². The heavy result
is tail-limited here. Direct quadrature samples are separate
diagnostics, not inputs in the proof.

## Keep the visible pair source, heat and reverse loss

Now allow the [vacuum charged-lepton inputs and elastic heating](INDEPENDENT-WARM-LATE-LEPTON-RECYCLING-2026-11-26.md)
and also dark fermion annihilation to their charged-lepton
pairs. Split the latter by a chosen invariant energy
`E_cut < m_V`. Below it, omit final Pauli blocking and use
the subpole denominator lower bound. The sum of the three
lepton species has `σv ≤ 3 C E_cut²/[4π(m_V²−E_cut²)²]`,
with squared vector coupling product C. Add this positive
coefficient to the internal quadratic loss ceiling.
The lepton final mass factor is at most one, and the
incoming mass factor is at most 3/2.

Above that invariant cut, each event removes two fermions
but extracts at least `E_cut` of dark energy. Nonnegative
pressure, the prior pair-source energy upper and the
Gronwall heat upper bound **gross** energy available for
extraction by `E_avail`. This remains true if extraction
and positive injection occur together: integrate the
energy inequality with the negative extraction term kept.
Thus the high-energy loss integrated in comoving units
is at most `Z_max=2E_avail/E_cut`, independently of how
large a near-resonant rate may be. Adding the actual
cumulative high loss Z back to X produces a quantity
`X+Z` whose derivative is at least the pure negative
quadratic evaluated at `X+Z`. Subtract `Z_max` from
that simpler solution. Each chosen cut alone gives a
lower bound; take the strongest of the tested cuts.

| q32, uniform F_B | Isolated fermion Ωh² lower | With all included vacuum charged-lepton terms |
| --- | ---: | ---: |
| 0 | 25.6926 | 24.6465 |
| 1 | 8.3270 | 7.3510 |
| 2 | 3.9158 | 2.9633 |
| 4 | 1.4529 | .5171 |

Under the vacuum lepton comparison, the imposed cap
would have to allow at least 4.89 to stop excluding a
total .120 by this bound; for the isolated case the
corresponding number is 16.54. A weaker imposed cap does
not demonstrate successful depletion. No dynamical scalar
occupation ceiling follows, and a bound on the outgoing
scalar occupation from a ballistic branch must not be
imported into this more general problem. Hadronic loss,
other interactions, and changes to the visible history
are not bounded here.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_fermion_depletion_kernel_reference.py > /tmp/fermion_kernel.json
python -W error warm_fermion_survival_with_leptons_reference.py /tmp/fermion_kernel.json \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-late-visible-recycling-envelope-2026-11-26.json
~~~
