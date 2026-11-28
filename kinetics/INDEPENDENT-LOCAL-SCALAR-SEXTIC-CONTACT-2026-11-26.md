# Low-momentum scalar sextic contact as a number-source control

Use the illustrative zero-momentum one-loop coefficient
c₆≈2.12e-16 GeV⁻² of the local interaction c₆φ⁶/6!, as matched in the
[upstream scalar completion](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-SCALAR-QUANTUM-COMPLETION-2026-11-24.md).
It is not the exact relativistic Yukawa-loop amplitude, and a different
UV matching condition can supply a different local coefficient. Use
5 or 10 MeV maximum center-of-mass energy to illustrate the restricted
low-invariant domain. Contributions above these cuts from the full
loop are not bounded by this EFT calculation.

In vacuum the contact matrix element squared is c₆². Final identical
particle normalization 1/4! and the massless four-body upper phase
space Φ₄≤s²/[24(4π)⁵] imply
σv_CM ≤ c₆² s/[576(4π)⁵]. The lab Møller ratio is s/(4E₁E₂).
For isotropic incident spectra, integrate this upper expression
analytically in the initial relative angle up to both the physical
four-body threshold and the selected CM cut. A separate numerical
angle integration checks that polynomial integral.

At a 10 MeV CM cut, the unoccupied-final number gain relative to the
existing scalar inventory per Hubble time is bounded by 4.60e-35,
9.51e-35 and 3.27e-41 at the archived 12 MeV, 9.3 MeV and low-reheat
3.157 MeV snapshots, respectively. These are rates for the restricted
local term, not evidence that hard pairs have no effect.

An even more optimistic budget lets every retained dark particle be
a scalar, lets every pair use the maximal sub-cutoff rate, drops inverse
reactions, and keeps adiabatic visible entropy with a photon Hubble
floor. The 12 MeV archive then gives optical depth <1.63e-31
from the 10 MeV contact. Accounting for maximally constructive
interference with the [previous quartic upper bound](INDEPENDENT-SCALAR-QUARTIC-MULTIPLICATION-BOUND-2026-11-26.md)
does not materially increase that quartic ceiling. Within this toy
comparison, a necessary effective |c₆| to permit 99% of the observed
cold-dark-matter rest-mass budget is ≳0.52 GeV⁻², over 2e15 times
the illustrative loop matching value. A coefficient at that threshold
is not known to work, and the model outside the cutoff is untested.

Do not apply this calculation to a coherent field, arbitrary final
Bose populations, or to hard finite-momentum heavy-fermion loops.

```sh
OPENBLAS_NUM_THREADS=1 python -W error scalar_local_sextic_contact_gate_reference.py EOS.late.dat results/conditional-relic-normalization-2026-11-26.json results/{12-quantum-48-shape,early-9p3-quantum-shape,low-reheat-early-pair-vacuum-shape}-2026-11-26.json
```
