# One internally produced exterior daughter cannot reach the high heavy band

The [persistent-support comparison](INDEPENDENT-WARM-HEAVY-HIGH-BAND-FINITE-DAUGHTER-SUPPORT-2026-11-26.md)
excluded direct heavy creation above 1 GeV comoving momentum
from two retained daughters. Relax this support restriction
for a single internally generated exterior light fermion or
scalar. It may be produced from two retained core inputs
and interact once more with another retained core input.
It need not stay below the retained edge after its creation.

Let (bge1) be the scale factor relative to the anchor,
(P) the core comoving cut, (m=20) MeV the larger of
the two light rest masses, and (M=100) MeV the heavy
mass. Each core input contributes at most
(sqrt{P^2+b^2m^2}) to comoving energy. If a propagating
light/scalar intermediate with mass (m_o) has (k)
core input leaves (counted with multiplicity), its
creation limits its squared comoving momentum to
(k^2P^2+(k^2m^2-m_o^2)b_{m birth}^2).
The factor of (b_{m birth}^2) is nonnegative since
(m_ole m), so the same inequality using later (b)
continues to hold during collisionless redshift.
Energy conservation at every subsequent interaction,
dropping all other positive-energy products, adds the
incoming allowances. This inductively gives the
ceiling at heavy creation
[
q_H^2le k^2P^2+(k^2m^2-M^2)b^2.
]
If the right side is negative, that branch cannot
create a heavy state; the numerical tables display zero
as the clipped momentum ceiling.

For one exterior daughter produced from two core
particles and incident beside one more core particle,
(k=3). Since (3m<M), the ceiling cannot increase
at later expansion and its largest value is .907732 GeV
(q32) or .448547 GeV (q16). This branch supplies no
heavy creation above the 1 GeV cut despite giving
one input unretained support. More generally, q16 is
excluded with at most five such core leaves for all
(bge1), or six up through (b<6.20891). At q32
the bound does not exclude four leaves, including
two one-hop exterior inputs; failing to exclude them
is not a positive rate or proof of attainability.

The bound requires light/scalar intermediates, core
roots subject to the stated cutoff at their use, and
no hidden heating. It does not cover inherited cohorts,
visible injection, arbitrary repeated re-entry, or prior
heavy propagation or elastic redistribution. That missing
population is separate from the
[gross reacting-exterior budget](INDEPENDENT-WARM-HEAVY-PERIPHERAL-GROSS-BUDGET-2026-11-26.md).
This test is about where heavy particles can be created,
not about how much surviving light or heavy matter remains.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_heavy_peripheral_ancestry_reference.py \
 results/warm-heavy-high-band-finite-daughter-support-2026-11-26.json
~~~

[Diagnostics](results/warm-heavy-peripheral-ancestry-2026-11-26.json).
