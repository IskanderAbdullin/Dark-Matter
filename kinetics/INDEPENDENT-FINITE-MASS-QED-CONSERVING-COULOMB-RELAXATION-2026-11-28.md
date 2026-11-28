# Prescribed Coulomb relaxation with a conserved density

The [collisionless Coulomb test](INDEPENDENT-FINITE-MASS-QED-LONGITUDINAL-SUM-2026-11-28.md) provides a frozen response with an undamped timelike pole. To test how an external relaxation clock can broaden it while preserving charge continuity, use a number-projecting relaxation model in the sense of [Mermin](https://doi.org/10.1103/PhysRevB.1.2362). The clock ν here is prescribed, not fitted to photon or charged-particle collisions. Only the **retarded** longitudinal density response is computed. A charge-conserving model need not conserve its momentum or energy, and the calculation does not infer an interacting fluctuation law or a physical damping rate.

Let χ₀(ω,k) denote the collisionless density polarization with χ₀(0,k)=−m_D² in this sign convention. For a frozen isotropic response, normalize the streaming term by k and set z=(ω+iν)/k. The angular occupation perturbation at velocity v_z contains the drive proportional to v_z and the relaxation target proportional to i(ν/k) δμ. Its current moment and density moment use the same underlying jump weight. Choosing δμ so the relaxation term has zero density moment gives

```text
Πν(ω,k) = (ω+iν) Π0(ω+iν,k) /
          [ω + iν Π0(ω+iν,k)/Π0(0,k)].
```

An independent angular solution of this projected kinetic equation reproduces Πν to 2.38e−14 and the continuity identity k δj = ω δn to 3.80e−16. Projecting against conserved total charge versus separately projecting each of the two charged species gives distinct results; the checked three examples differ at most 2.37e−8 relatively. They are alternative assumptions, not proven equivalent. The displayed pole and sum results use the total-charge projection and equal ν for the two species.

The static Debye value is fixed at any nonzero ν by the conserving denominator. A simple substitution ω→ω+iν without that denominator changes the static polarization, by as much as 55.4% in the sampled controls. The unsubtracted Coulomb response still has the instantaneous high-frequency piece 1/k²; the tested dispersive moment is its difference from the static response. Both this moment and three imaginary-axis dispersive values are checked against direct response. The first frequency-weighted moment is independently checked against the collisionless leading high-frequency coefficient inside the same kinetic model. This last agreement is not a microscopic many-body sum-rule proof.

| Check | Worst defect |
| --- | ---: |
| Static inverse-frequency spectral weight | 6.54e−13 absolute |
| Imaginary-axis normalized spectral reconstruction | 4.38e−10 absolute |
| Frequency-weighted moment, divided by the model coefficient | 4.25e−9 |
| Small-k complex pole vs the Drude limit | 1.09e−9 relative |

For each archive the three input ν values are 1e−5, 1e−4 and 4e−4. The 30 tested root continuations stay on the shifted analytic sheet: 0 < −Im ω < ν. At k=1e−7 the ratio −Im ω/ν is about .5, as in the small-k Drude limit. For the ν=1e−4 fine test it grows to .530 at k=.0008, .753 at .00256, and .931 at .004. This sensitivity to phase speed is not a calculation of collisional production. No extension to a transverse relaxation kernel, physical plasmon population, scalar source, or matched rate is made.

Reproduce with `python -W error finite_mass_qed_conserving_coulomb_relaxation_reference.py` from `kinetics/` and compare the [receipt](results/finite-mass-qed-conserving-coulomb-relaxation-2026-11-28.json).
