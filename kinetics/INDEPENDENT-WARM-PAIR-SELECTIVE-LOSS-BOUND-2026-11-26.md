# Add a momentum-moment constraint to the late support bound

This is a **narrower** comparison to the existing
[finite-support bound](INDEPENDENT-WARM-PAIR-LATE-SUPPORT-BOUND-2026-11-26.md).
Start at the archived 1.206 MeV state. Allow the forward pair reaction
to remove any retained parent, but assume the survivors only redshift:
no inverse gains, elastic redistribution, new incident population, or
re-entry. The previously excluded count is inert. Apply the same
assumed photon floor and a cap F on future final-state Bose occupation.

The retained comoving sum of squared momenta can therefore only
decrease, even if the particles are lost selectively. The moment per
*survivor* need not decrease, and is not assumed to. Bound its
initial sum with the exact shell kinetic-energy moment, using
(p^2 = (2m+K)K) and the shell extrema of (K).
The mean p² is in ([5.1033,5.3314]\times10^{-5}) GeV²;
direct positive quadrature gives a value inside that interval.

For uncorrelated isotropic incoming directions,
`<p_cm²>_angle = (E1 E2 - m²)/2 <= (p1² + p2²)/4`.
Combine this with the pointwise finite-support estimate and the
same amplitude/flux inequalities as the previous bound.
At equal charge populations let x be the surviving fraction of
the **retained** stable count, and let r be the initial mean p²
upper bound divided by 2 P². In a time variable μ normalized to
the integrated support bound, the largest allowed loss obeys
`-dx/dμ <= min(x², r x)`. With arbitrary isotropic charge
fractions the conservative substitution is twice r. Until
(x=r) the comparison decays exponentially; below that point it
follows the remaining support-only inverse-count solution. This
combination is stronger than applying either bound alone.
A separate integration of the comparison ODE checks both regimes.

At the nominal photon floor and equal charges:

| Assumed future occupation cap F | Conditional lower stable Ωh² |
| ---: | ---: |
| 0 | 56.06 |
| 0.001 | 56.06 |
| 1 | 52.10 |
| 2 | 46.11 |

At F=1 the analogous support-only lower value is about 27.12.
In this narrower comparison a cap must allow F of at least 43.66
before the lower envelope reaches 0.120. With arbitrary charge
fractions and a 2% lower photon temperature floor it must still
allow at least 42.28. These are necessary conditions to reach
the target within the assumptions, not estimates of the actual
occupation, and not sufficient to achieve that count.

This moment step **cannot** be justified by an anchor snapshot
alone if inverse reactions or scattering can repopulate momenta.
It does not supersede the less restrictive finite-support bound,
or exclude the full model or any other production history.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error warm_pair_selective_loss_bound_reference.py \
  results/warm-coupled-pair-b10-28-shape-2026-11-26.json \
  results/warm-coupled-pair-b10-count-audit-2026-11-26.json
~~~
