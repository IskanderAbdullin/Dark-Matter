# Close the Bose cap in a strictly loss-only warm continuation

The [late pair bounds](INDEPENDENT-WARM-PAIR-ENERGY-AND-BAND-2026-11-26.md)
take a future scalar occupation cap as an assumption.
Here it follows from a comparison calculation in a **narrower** problem.
Start from the reconstructed 1.206 MeV retained shape. Parents
can redshift and be destroyed by forward pairs only. There are
no inverse gains, redistributing collisions, other sources,
or initial exterior scalars supplying occupation. The old
stable exterior count remains inert. The newly produced scalar
population need not remain inside the archived mesh. No cap
on its later momentum support is used in the gain estimate.

The maximum per-mode parent occupation cannot increase along
the stable characteristics; its archive value is (S=5.725\times10^{-6}).
Incident stable momentum is below (P/b).
The scalar lab energy at any achievable final state exceeds
(E_{\min}=\sqrt{m^2+P^2}-P=9.087\) MeV, or a greater value later.
The same amplitude bound \(C(P)\,p_{\rm cm}^2\) applies.

For a fixed outgoing scalar momentum, use the spatial delta
function to eliminate its partner's momentum. Integrate the
remaining energy delta over the second incident momentum's
polar angle around the difference of the first incident
momentum and the fixed outgoing momentum. The angular
Jacobian cancels the partner-energy denominator and leaves no
worse singularity than the inverse norm of that difference.
For a ball of incident radius (P/b), its 3D integral is
bounded by (2\pi(P/b)^2), including when the fixed momentum
is inside the ball. Bound the remaining parent radial integral
by ((P/b)^2/(2m)). Thus the rate for any outgoing occupation
is no larger than

~~~text
A(b) (1+F)^2,   A(b) = G C(P) P^6 S^2 / [32 E_min (2*pi)^3 m^2 b^6],
~~~

where F is the global scalar occupation supremum and G
overcounts spin/label factors. Using G=32 deliberately
overcounts the product of the model's stable degeneracies
and even includes both possible final labels without the
identical-pair reduction. The photon floor implies
`integral A dt <= A(1)/(4 H_gamma)`.
Comparison along every scalar characteristic gives
`1/(1+F) >= 1/(1+F_initial) - integral A dt`,
without differentiating the supremum.

With (F_{\rm initial}=7.0562\times10^{-4}), G=32 and a
2% lower photon-temperature floor, the resulting global cap
is **0.002109** (rounded up). The stronger loss-only
moment-plus-support comparison then gives a stable fixed-count
Ωh² floor of about **56.00**, not 0.120. At the nominal floor
the cap is 0.002053 and the abundance floor is 56.055.
The alternative spin multiplier and photon floor are listed
separately in the archived results.

This does **not** close the Bose cap for the full signed model.
Inverse events or elastic scattering can replenish the parent
occupation; unspecified older exterior scalars could also
stimulate the process. The calculation is a conditional
failure of forward-only ballistic depletion to reach the
target, not a late-time prediction for a coupled dark sector.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error warm_pair_bose_gain_ceiling_reference.py \
  results/warm-coupled-pair-b10-28-shape-2026-11-26.json \
  results/warm-coupled-pair-b10-count-audit-2026-11-26.json
~~~
