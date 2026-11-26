# Count budget for a dominant 86 eV component

This comparison is independent of the late pair rate, stimulation
cap, momentum grid and transfer solver. Use the already archived
rest-density normalizations for the two physical masses: a
20 MeV stable fermion and an 86 eV scalar. Suppose that after
the stated checkpoint the **total comoving number in these
two species** and the visible comoving entropy are conserved.
Equal-number conversion (\chi\bar\chi\leftrightarrow\phi\phi)
and elastic redistribution satisfy the number condition.
Additional production, number-changing reactions, transfer
to other species, or a different entropy history need not.

Write (r=m_\phi/m_\chi) and
`Y = Omega_phi h² + r Omega_chi h²`. This is the sum of the two
counts expressed as scalar rest density. It is unchanged by
closed conversion, independent of how its efficiency varies
across momentum or time. Consequently:

- Converting every available stable particle gives no more
  scalar rest density than Y.
- If the total rest density is fixed to the reference .120,
  the only possible closed partition has
  `Omega_phi h² = (Y - r * .120) / (1-r)`, when feasible.
  This equation specifies a count budget, not a rate or a
  guarantee the required partition can be attained.
- For half the .120 target to be scalar, the available number
  per visible entropy would have to give
  `Y >= .060 * (1+r)`.

| Archived history | Maximum scalar rest density from its count Y | Scalar share if total rest density is .120 | Required count/entropy increase for a 50% scalar share |
| --- | ---: | ---: | ---: |
| Low reheat, instantaneous | 5.1631×10⁻⁷ | 2.56×10⁻⁹ | 1.162×10⁵ |
| Low reheat, finite Tmax/TR = 5 | 5.1875×10⁻⁷ | 2.29×10⁻⁸ | 1.157×10⁵ |
| Warm 50 MeV branch, q16/q32 | 5.2500×10⁻⁴ | 0.4371% | 114.3 |

The low-reheat inputs are the
[shifted prompt spectrum](INDEPENDENT-FINITE-COOLING-SHIFTED-COUNT-2026-11-26.md)
and its early forward addition at the same 1 MeV anchor. The
added scalar number and removed stable number cancel exactly in
Y. In the [warm q16/q32 source](INDEPENDENT-WARM-SIGNED-Q32-PRODUCTION-2026-11-26.md)
the old retained and outside counts are used together in Y,
without inventing old exterior momentum shapes. For the
[evolved q16 comparison to 1.206 MeV](INDEPENDENT-WARM-COUPLED-B10-2026-11-26.md)
the reconstructed Y agrees with its own 12 MeV budget within
6.3×10⁻¹⁴ relative, even though the separate stable count changes.

This closes a specific loophole: faster **equal-number** conversion
of any of these archived populations cannot by itself account
for a dominant 86 eV rest component at the specified total.
Increasing the new scalar count or changing the normalization
requires leaving the closed inventory assumptions. The result
does not estimate extra production, bound all other mechanisms,
or establish an observational exclusion.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_and_low_closed_number_budget_reference.py \
  "$D"/finite-cooling-shifted-dist-2026-11-26.json "$D"/finite-cooling-shifted-spectrum-shape-2026-11-26.json \
  "$D"/50MeV-from-empty-pair-quantum-q32-audit-2026-11-26.json "$D"/warm-coupled-pair-b10-count-audit-2026-11-26.json
~~~
