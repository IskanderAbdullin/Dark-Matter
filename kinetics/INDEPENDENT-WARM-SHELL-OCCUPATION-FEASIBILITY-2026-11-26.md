# What saved warm fermion shell moments can tell us about an occupation cap

At the six matched saved times from about 10.391 to 10.137 MeV,
the measured peak of the retained light fermion **reconstruction**
is of order 6e-5 per spin. Test a different question: without
assuming the entropy reconstruction, can the archived light
number and kinetic energy of each shell force or rule out
occupation above 1/2? This matters because the
[signed survival comparison](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
uses a pointwise cap valid at all times.

Use radial measure (dV=p^2 dp) and monotonically increasing
(K(p)=sqrt{p^2+m^2}-m). For cap (c), density is feasible if
(0le Nle cV) and the kinetic-energy moment lies between the
low- and high-edge packed profiles. Each packed profile fills
exactly volume (N/c) at height (c). Exchanging occupancy from
higher to lower (K) can only decrease the moment. Convexly mixing
the two packed profiles produces all intervening energies at the
same number and respects the cap. Thus these conditions are both
necessary and sufficient, up to numerical quadrature in this audit.
Bisection brackets the smallest feasible cap, which is a **lower**
bound on the true essential maximum, not an upper bound.

| Momentum order | Range of smallest compatible `max f` across the six endpoints | Range of peaks of the actual entropy reconstructions |
| --- | ---: | ---: |
| q16 | 5.91477e-5 to 6.02327e-5 | 5.91578e-5 to 6.04338e-5 |
| q32 | 5.77992e-5 to 5.84644e-5 | 5.78493e-5 to 5.84890e-5 |

All tested shell moments admit low-occupation profiles below 1/2.
But at **every** saved endpoint, the script also constructs a
shell with (f=1) on a narrow subinterval and zero outside it,
with the same two moments in that shell; leave all other shells
unchanged. Its occupied volume sets the number exactly. Sliding
that interval from left to right changes its energy continuously
from the cap-1 lower to its upper extremum. Independent quadrature
checks reported relative errors below (7	imes10^{-14}) for
the recovered input moments and (3	imes10^{-16}) for the
witness kinetic energies. A mathematical witness is not evidence
that such a sharp profile forms dynamically. It shows that
these two saved moments alone cannot establish the pointwise cap.

The massless negative control ([0,1],N=.1,K_{m moment}=.09)
passes the cap-1 test but exceeds the cap-1/2 energy limit
.08815969. A near-edge dilute control avoids subtracting
indistinguishable floating-point cube roots; it rationalizes
the width of the occupied interval before integrating.

This audit says nothing about unresolved exterior cohorts, the
heavy state's occupation, other times, or the actual continuous
collision dynamics. A physical cap would require evolution of
the local distribution, or a separate regularity or comparison
bound not furnished by the two moments.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_shell_occupation_feasibility_reference.py
~~~

[Machine-readable diagnostics](results/warm-shell-occupation-feasibility-2026-11-26.json).
