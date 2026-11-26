# Update the conditional fermion survival comparison at 9.748 MeV

The [selected-pion signed comparison](INDEPENDENT-WARM-SIGNED-FERMION-WITH-CHARGED-PIONS-2026-11-26.md)
used a single global radiation-clock lower floor after its
10.442 MeV anchor. Its key hypothesis is that **both** light and
heavy fermion occupation per spin remains below 1/2. The
scalar occupation is not capped.

Preserve the same signed reaction and selected-visible
bounds. On each tabulated radiation-EOS interval take
the independent worst-case endpoints of the entropy and
energy degrees of freedom (the interpolant is monotone on
each interval). Split the integration grid at those knots.
Integrate each expansion power exactly against its local
lower Hubble floor, with a photon-only late-time fallback.
Every integrated upper rate weight is independently checked
against the original global-floor weight. This changes
the original conditional ( Omega h^2) lower result from
1.119225 to 1.121197 at the fine resolution, a small correction.

Separately, condition on the measured [q64 checkpoint at
9.748 MeV](INDEPENDENT-WARM-Q64-HISTORY-SERIAL-2026-11-26.md)
as the new starting state. Its represented fermion count
is the lower initial count; do not credit an unknown heavy
or old exterior number. The original full future available
energy ceiling already charges the selected injection and
heating budgets. Transport that ceiling with (b^3ho)
from the older entropy anchor (relative (b=1.071408));
do not add the same old budget twice. The number-to-density
normalization uses the same archived visible entropy
convention and physical conversion. The revised conditional
lower result is 2.111328 at cap 1/2, versus the target .120.
With all rates otherwise fixed, the later comparator still
certifies exclusion for a larger available energy ceiling,
until roughly 1.267 times the specified budget. This last
number measures a sensitivity in the comparison, not a
bound on unknown new physical sources.

Neither the actual half-occupation barrier nor the validity
of the later starting spectrum in the fully coupled physical
history is proved. The numerical q64 reconstruction alone
cannot establish a mode cap, a complete hadronic network, or
a cosmological exclusion without the stated assumptions.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_signed_fermion_piecewise_clock_reference.py EOS.dat
~~~

[Diagnostics](results/warm-signed-fermion-piecewise-clock-2026-11-26.json).
