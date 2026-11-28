# How much recorded exterior birth energy can be above the high-heavy input threshold?

The [gross exterior requirement](INDEPENDENT-WARM-HEAVY-PERIPHERAL-GROSS-BUDGET-2026-11-26.md)
needs the energy of exterior inputs that actually react, not a
one-way receipt. As a narrower diagnostic, identify an energetic subset
of the stored receipts that could supply high heavy creation, even when
**both** inputs are exterior.

For a heavy creation cut (Q=1) GeV comoving, a light or scalar input
with (q<Q/2) and mass below (M/2) always carries less than half the
energy of any high heavy output. A two-input event with both below that
threshold cannot create a high heavy output. If one input is below the
threshold, the other covers more than half the energy of the high
heavy outputs summed over the event; if both are above, their total
covers all of it. Thus a ledger of **gross** energy from inputs with
(qge Q/2) would suffice with a factor of two for this part of the
high-band production budget. This leaves visible inputs and heavy
redistribution separate.

A stored high-side cohort has number (N), mean birth energy (ar E)
and guaranteed lower birth energy (E_b). For birth threshold (E_*>E_b),
the count above threshold is at most (N(ar E-E_b)/(E_*-E_b)), capped
by (N). Because (E/(E-E_b)) decreases with (E), the energy above
threshold is at most (N E_*(ar E-E_b)/(E_*-E_b)), capped by (Nar E).
A floor-plus-threshold two-point control attains both bounds in the
nontrivial regime. Use each cohort's own recorded birth expansion to
convert (q=Q/2) to (E_*), without treating the mean as an individual
momentum cutoff.

From the 10.442 MeV restart to (b=1.0302) the tracked **new** q32
cohorts have at most (1.98	imes10^{-24}) GeV⁴ of eligible birth
energy (about 7.61% of their tracked outer birth energy), in anchor
comoving number times physical birth energy units. The corresponding
q16 upper is (2.17	imes10^{-18}) GeV⁴ (about 8.19%). Stored inherited
q32/q16 receipts add upper values of (1.06	imes10^{-23}) and
(1.11	imes10^{-17}) GeV⁴, respectively. The audited receipts contain only upper-side cohorts. A lower-side cohort stays
below the tag if its support edge is below the tag and there is no reheating.

This is an upper on an eligible **one-use** energy pool under
collisionless propagation before use. It is not the true reacting
energy, and it misses old untracked cohorts, acceleration of
initially low-momentum input, repeated use and later new production.
In particular it is not a complete closure of the high heavy band.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_heavy_eligible_exterior_receipt_reference.py
~~~

[Diagnostics](results/warm-heavy-eligible-exterior-receipt-2026-11-26.json).
