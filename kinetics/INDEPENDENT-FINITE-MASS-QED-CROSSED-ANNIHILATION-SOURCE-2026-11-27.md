# Hard crossed annihilation is a small sink in a frozen cold-mode comparator

The [crossed-annihilation threshold](INDEPENDENT-FINITE-MASS-QED-CROSSED-ANNIHILATION-2026-11-27.md)
excludes this channel from the cold mode's soft photon window, not from its
full collision operator. Integrate **one** charge assignment
(parent + anti-daughter ↔ scalar + photon) at scalar momentum 4.1e−6
using the two frozen archive spectra. Use the daughter histogram as an
explicit **proxy** for the anti-daughter occupation, the held scalar
occupation, and a photon Bose distribution at temperature 1 in the
archive's dimensionless units. The inverse term uses both Pauli factors.

Use unaveraged spin and polarization sums, without a second charge
assignment or chemical equilibration assumptions. After integrating the
photon's azimuth in an isotropic bath and averaging the internal pair
azimuth, the measure at fixed scalar momentum is
$k r_* /(128\pi^3 W E_s)$ per photon energy and the two polar cosines,
where $W$ is the pair invariant mass and $r_*$ its COM momentum.
A separate unit-amplitude inverse cross-section times Møller flux
reproduces the phase normalization. The Born comparator uses the same
one-charge conventions and the external scalar collision clock,
including an identical Yukawa coupling squared which cancels in the
ratio. Only the crossed channel carries the additional photon coupling
squared.

The exact crossed spin sum decomposes into two positive-denominator
first- and second-pole terms in the two charged-photon invariants.
Integrate their azimuth and then each constant occupation interval in
the boost polar cosine analytically. Independent angular quadrature and
an explicit four-vector boost validate this reduction. Radial
quadrature resolves the pair threshold, transformed histogram edges,
and crossings where two occupation boundaries exchange order. The
finite histogram gives no forward source above the sum of the two
maximum charged energies minus the fixed scalar energy. A phase-space
bound on the remaining inverse Bose tail beyond photon energy 68 is
below 6.8e−24 in the unscaled loss integral in both controls. With
FD/BE weights all evaluated at a common temperature, each of the
four integrated tests obeys detailed balance.

| Archive | Crossed net / Born net at g=.005 | Anti-daughter shape scale required for zero | Pauli-allowed maximum shape scale |
|---|---:|---:|---:|
| fine | -2.20520e-06 | 11.6832 | 9.3840 |
| third | -2.53590e-06 | 11.3478 | 9.5764 |

Final mesh drift across gain, loss, and both net controls is at most 4.83e-12 relative. The angular control is within 3.11e-12; equilibrium net/turnover is at most 1.63e-16. The proxy therefore remains a sink under every Pauli-allowed uniform rescaling of the stored anti-daughter shape.

These are *frozen single-mode* corrections, not integrated abundance,
nor an inferred antiparticle distribution. In particular the fitted
neutrality scale only multiplies the stored daughter shape; it does
not exclude other anti-daughter distributions. Another archive time,
momentum, charge state, or photon temperature can change the ratio
or sign. No thermal vertex, alternative hard channel, full reaction
network, or cosmology is inferred.

Replay in kinetics with the two pinned spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_crossed_annihilation_source_reference.py
~~~
[Receipt](results/finite-mass-qed-crossed-annihilation-source-2026-11-27.json).
