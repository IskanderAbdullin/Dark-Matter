# The lowered halo produces entropy even in collisions that eject nobody

The [escape ledger](INDEPENDENT-HALO-BOSE-ESCAPE-ENERGY-LEDGER-2026-11-26.md)
tests direct loss across a selected energy cut.
This check removes those events and asks whether
the lowered distribution is a local fixed point
of the remaining elastic collisions. Include a
collision only if both incoming and both outgoing
energies lie below the same tidal threshold.
This subset is closed under reversing the collision and
preserves number and pair energy. It is a diagnostic
restriction, not a physical rule that forbids
escapes in a real halo.

For the quantum gain and loss weights
$G=f_3f_4(1+f_1)(1+f_2)$ and
$L=f_1f_2(1+f_3)(1+f_4)$, local production
has the nonnegative integrand
$(L-G)\ln(L/G)/4$, with the conventional total
elastic cross-section measure and both labeled
incoming momenta. For importance sampling using
the incoming $f_1f_2$ weights, evaluate it as
$L\Delta\tanh(\Delta/2)/2$, where
$\Delta=\ln(L/G)$. This uses microreversibility
and the symmetry of $\Delta\tanh(\Delta/2)$
under reversing a collision. An independently
sampled sign-changing form $L\Delta/2$ gives
the same integrated result. Direct use of
$(1-G/L)\Delta$ with only $L$-side sampling
had badly weighted rare events; the bounded
hyperbolic factor avoids that estimator defect.

An unlowered Bose thermal form restricted to the
same elastic allowed collision pairs would give
zero production: $\ln[f/(1+f)]$ is an affine
function of energy. Subtracting the finite value
of the population at escape makes it nonlinear,
so the lowered state is not that local equilibrium.

| Escape depth, peak occupation 7,000 | Initial internal entropy production per scalar per Gyr |
|---:|---:|
| 2 | .008149 |
| 5 | .002012 |
| 8 | .0001945 |

Entropy uses the dimensionless $S/k_B$ convention.
All examples retain the illustrative
$1.25\,{\rm cm^2/g}$ vacuum cross section.
The independent symmetrized form agrees to within
.02%; the combined grid/scramble effect is up to
1.06%. This does not count entropy carried by exported
particles, external heat, or evolution of the
density and velocity distribution.

~~~sh
python -W error halo_bose_lowered_uu_entropy_reference.py
~~~

[Initial internal entropy ledger](results/halo-bose-lowered-uu-entropy-2026-11-26.json).
