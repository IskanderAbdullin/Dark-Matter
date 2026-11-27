# The same late exterior receipts can give either sign of instantaneous high-tail flow

At the exact corrected fine q64 checkpoint (7.394049 MeV),
reactivate the 536 known collisionless receipts in three frozen
two-packet realizations. In each cohort put the low packet
between the original birth floor and mean, and set the second
to a prescribed physical momentum of .68, .80, or 1.00 GeV
at the late anchor. Preserve the count and mean **birth**
energy of every cohort. Some low packets now fall inside the
retained mesh and are added there, rather than excluded or
counted twice. Recover the new shell number and kinetic moments
using positive reconstructed occupations. Empty the old
collisionless ledger only in these synthetic snapshots.

The restored total energy is only 2.826e-11 of the old retained
energy. At a fixed physical energy cut C=.71482 GeV, the
recorded-cohort [ballistic energy ceiling](INDEPENDENT-WARM-CONVEX-HIGH-ENERGY-TAIL-2026-11-26.md)
is 1.437e-24 GeV⁴ above the cut. The point-packet realizations
carry respectively 0, 1.102e-24 and 1.026e-24 above it.
This is a check on the unprojected moments, not a claim that
any one is the true exterior shape.

| Hard packet at the anchor (GeV) | Gross inverse parent energy | Signed production of excess above C |
| ---: | ---: | ---: |
| .68 | 2.535e-6 | +1.478e-10 |
| .80 | 1.008e-5 | −9.489e-7 |
| 1.00 | 5.960e-6 | −1.222e-6 |

Rates are normalized to the represented dark energy times
H, for the open prompt angular diagnostic with stationary
parent blocking. With the .68 GeV choice, the daughters
initially lie below C but some outgoing daughters can cross it.
In the other two choices the already present high energy is
redistributed downward in the instantaneous signed balance.
Halving the hard shell half width from .005 to .0025 GeV
preserves all three signs. The .68 rate changes by about 6%;
the other two signed rates change by less than 0.03%.

A small stored energy alone therefore cannot be used as an
instantaneous rate ceiling without decay/return dynamics or
a phase constraint. These are deliberately specified frozen
shapes, not a dynamical upper on tail growth or its lifetime.
In particular, multiplying any displayed rate by a Hubble
time is not a justified exposure estimate.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for w in .005 .0025; do
  python -W error warm_corrected_exterior_reactivated_tail_stress_reference.py EOS.dat \
    results/warm-corrected-anchor-signed-floor-2026-11-26.json \
    results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz --half-width "$w"
done
~~~

[Base shell width](results/warm-late-reactivated-tail-width005-2026-11-26.json);
[half-width control](results/warm-late-reactivated-tail-width0025-2026-11-26.json).
