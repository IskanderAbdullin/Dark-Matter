# The charged bath creates a transverse photon gap and moves the absorption opening

The preceding [matched real support probe](INDEPENDENT-FINITE-MASS-QED-SIGNED-TRANSFER-MATCH-2026-11-27.md)
uses **vacuum lightlike photons**. At those energies a necessary
check is whether a collisionless transverse electromagnetic pole
exists. As a separate controlled model, freeze the two isotropic
charged archive distributions and use the leading weak-coupling
relativistic kinetic polarization with their nonzero masses.
Give each charged species two spins, equal absolute coupling g=.005,
and the stated particle-shaped antiparticle proxy. No collisional
width or photon occupation is assumed for this retarded response.

The transverse inverse propagator is
$\omega^2-k^2-\Pi_T$, where $\Pi_T$ weights the derivative of each
charged occupation by $v^2$ and the angular average
$\frac14\int_{-1}^1(1-u^2)/(1-k vu/\omega)\,du$.
At $k=0$ the angular factor is $1/3$. The frozen histogram derivative
can be integrated exactly as its jumps or by parts as a positive
radial integral. These forms agree within 9e−16 on the two fixtures.
An independent angular integral, the finite-mass equilibrium rest
formula, and the small-$k^2$ pole expansion provide separate controls.
The standard equilibrium comparison is [Braaten and Segel](https://arxiv.org/pdf/hep-ph/9302213);
using the *archived nonthermal shapes* here is a frozen
collisionless kinetic model, not their thermal state.

| Archive | Charged-only Debye mass | Transverse frequency at k=0 | Scalar energy held fixed |
|---|---:|---:|---:|
| fine | .00288296 | .00166386 | .00087207 |
| third | .00288534 | .00166524 | .00086725 |

The positive radial form makes the timelike response at least
its zero-momentum value; the transverse pole has no energies
below this gap. Thus the vacuum lightlike real spectrum below the
gap cannot be counted as propagating transverse quasiparticles.
This does not rule out spacelike spectral weight.

For absorption the relevant charged transfer is $q-K_T$.
At small photon momentum it has negative energy and is timelike.
A heavier future-directed parent minus a lighter future-directed
daughter cannot equal such a transfer. As the transverse momentum
increases, the transfer first becomes spacelike at momentum
about .00127–.00131, depending on angle and archive (transverse
energy about .00215–.00218). This is a **past-null** opening,
not the earlier [future-null vacuum example](INDEPENDENT-FINITE-MASS-QED-ABSORPTION-CROSSES-NULL-TRANSFER-2026-11-27.md):
the required parent energy diverges as the opening is approached
from above. For transverse cosine zero the finite archived parent
and daughter supports first become *kinematically possible* at
momenta .001404/.001433 (fine) and .001423/.001452 (third).
Those edges do not estimate the rate or whether the tiny tail
there contributes appreciably.

A sensitivity control removes only the mirrored antiparticle
contribution to the kinetic polarization. It changes the gap to
.001177–.001178 and reduces the opening momentum, but retains the
low-energy transverse prohibition in this specific model. This is
not an independently neutral charged-plasma solution.
No polarization-resolved production matrix, longitudinal contribution,
Landau cut rate, finite virtual remainder, or cosmological transfer
follows from this gate. Replacing the vacuum photon energy by the
pole energy is insufficient to obtain those quantities.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_frozen_transverse_plasmon_gate_reference.py
~~~
[Raw receipt](results/finite-mass-qed-frozen-transverse-plasmon-gate-2026-11-27.json).
