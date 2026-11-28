# Frozen hard-pair support budget with exact archive-cell splits

The [bounded live hard test](INDEPENDENT-EARLY-POLE-LIVE-JOINT-HARD-PAIR-PACKETS-2026-11-26.md)
evolves photons only up to q=2 and outgoing charged particles up to
q=2.2. A direct estimate of what it omits needs both number and
energy, using the earlier resolved charged archives rather than the
synthetic initial hard cells.

At frozen x=.05 take the archived, piecewise-constant parent and
daughter occupations on each of two grids and held T_ref Bose photons.
Use the same spin-weighted vacuum central cut (s/T_ref²≥1,
|cos(theta*)|≤.8 and incoming q≥.5). On this one reference, the
incident upper bound is 40. Split the incoming quadrature at every
reported window boundary. At fixed center-of-mass scattering angle
solve every azimuth crossing of either archive edge and the outgoing
support edges. Use the common split for both grids so mismatched
classification does not create a false convergence difference.
These are signed net gains, not the probability of an event.

| Evolved support | Parent count / heat share | Daughter count / heat share |
| --- | --- | --- |
| Photon upper 2, fermion upper 2.2 | 46.6–47.2% / 32.6–33.2% | 26.4% / 15.1% |
| Photon upper 4, fermion upper 5 | 89.7–90.1% / 82.9–83.6% | 74.1% / 63.9% |
| Photon upper 8, fermion upper 10 | 97.9–98.0% / 98.6–98.7% | 91.0% / 92.2% |

All these evolved-support comparisons keep the outgoing lower
boundary q=.45. At least one daughter goes below .45 in a subset
carrying about 8.3% of the signed daughter count and 6.1% of the
signed daughter-channel heat. Almost all remaining daughter
omission at the widest upper support is therefore on the low
side. Parent low-side shares are about 1.9% and .9%.
Heat is the signed half-total event energy in the charged-count
normalization, with positive sign into the charges.

The last two split refinements change the compared shares by at
most about .04 percentage point. The old method using a sharp
indicator on an unsplit global quadrature did not converge for
these fractions and should not be substituted. At the strongest
illustrative gauge coupling in this check, the full retained vacuum
subset with held unit-temperature photons gives initial charged
count gains about 32–33 (parent) and 171 (daughter) in the
reference per-log-scale normalization. This is a frozen bath
comparison, not a calibrated cosmological rate: no evolving
photon bath, true medium amplitude, remaining soft phase space,
or sustained production history is included. A useful next
question is whether moving the *live* support toward the second
or third row stabilizes the evolved moment response.

Reproduce with early_pole_archive_hard_pair_exact_support_reference.py
(--level 0, 1, or 2). Raw:
results/early-pole-archive-hard-pair-exact-support-2026-11-26.json.gz.
