# Photon heat and scalar number require different real-channel moments

Reusing the full-support [charged-pair crossing](INDEPENDENT-FINITE-MASS-QED-CROSSED-ANNIHILATION-SOURCE-2026-11-27.md)
and [parent plus scalar crossing](INDEPENDENT-FINITE-MASS-QED-REAL-SCALAR-ABSORPTION-2026-11-27.md),
compute photon energy, scalar energy, and charged energy separately at
the same frozen scalar momentum 4.1e−6. Keep one charge assignment.
For the charged-pair process, use the same explicit daughter-shaped
anti-daughter proxy. No antiparticle proxy enters the second process.
The photon bath is Bose-distributed at temperature 1 in archive units.

In the scalar-producing direction the first reaction produces a
photon, while the second consumes one. Weight forward and inverse
histories with their *own* photon energies and signs. The scalar
work is $E_s$ times the scalar number source; the charged work
uses the event charged energy gain, not $E_s$ times a proxy number.
The signed energy ledger closes to the floating-point accuracy
reported in the receipt. A common-temperature scalar absorption
fixture gives no photon heating.

| Archive | Both channels: photon heat / Born scalar energy source | Both channels: scalar number / Born scalar number |
|---|---:|---:|
| fine | -0.00735729 | -2.67024e-06 |
| third | -0.00851323 | -3.081e-06 |

Both ratios use g=.005. Maximum two-mesh heat drift 3.81e-09; scaled energy ledger defect 1.26e-14; common-temperature heat/turnover defect 3.03e-16.

These are comparisons to the magnitude of a positive **Born scalar
energy source at this one mode**, not percentages of the total photon
bath energy or its full transfer rate. Both reactions decrease
the scalar source in these archived fixtures, but their photon
energy effects oppose each other. The net photon effect is much
larger in this normalization than the scalar number effect because
typical photon energies exceed the held scalar energy substantially.
No energy evolution, finite virtual correction or dressed heat
transport follows from this fixed-background two-channel control.

Replay with the pinned archived spectra from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_hard_real_heat_reference.py
~~~
[Fine receipt](results/finite-mass-qed-hard-real-heat-fine-2026-11-27.json),
[third receipt](results/finite-mass-qed-hard-real-heat-third-2026-11-27.json).
