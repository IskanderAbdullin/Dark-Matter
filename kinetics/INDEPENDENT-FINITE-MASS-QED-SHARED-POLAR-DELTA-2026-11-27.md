# One exact polar-delta kernel for the two positive charged cuts

The [scalar-latest contour comparison](INDEPENDENT-FINITE-MASS-QED-SINGLE-CUT-CONTOUR-2026-11-27.md)
found equal and opposite **leading** absorptive log coefficients
per occupation on two mapped positive charged cuts. In fact the
whole polar-delta kernel is shared at nonzero radial displacement,
not only its leading logarithm.

On simultaneous charged poles, the photon difference K obeys
Q_scalar·K=0. In the scalar rest frame the internal and external
charged momenta have the same magnitude p_star, itself common
to the parent and daughter. The remaining scalar-rest photon
virtuality is −2p_star²(1−cos(gamma)), with gamma the angle between
the two charged vectors. The azimuth integral has discriminant
(c_ext−c_int)², where the c's are their cosines relative to the
bath-frame boost axis. Since
delta_E=(|Q_scalar|/m_scalar) p_star (c_int−c_ext), its photon inverse
average is exactly −|Q_scalar|/(2 H |delta_E|),
H=m_scalar p_star. Combining the cut Jacobian and numerator
gives, with one overall statistical factor still omitted,

    common_kernel(delta_E) = −pi C/(4 |Q_scalar|)
                             + pi (P·b)/(2 H |delta_E|).

The remaining radial support is the common vacuum Born aperture.
Both lab angle derivations recover this kernel to 8.3e−15 relative
or better in the tested controls. Direct scalar-rest azimuth quadrature
agrees to 7.3e−15. Both polar roots leave the aperture on the same
outside tests.

Use the opposite contour signs to multiply the common kernel by
f_parent(E_internal)−f_daughter(E_internal−E_scalar). Under identical
smooth test occupations the independently combined cuts cancel
through the whole aperture, with residual below 1e−17 at the smallest
cutoff. For the unequal frozen archives, every cell integrates by a
linear term and a log. The competing occupations leave logarithm
coefficients .320569 (fine) and .242599 (third) at untuned centers;
the two deliberately tuned step events give .231915 and .142356.
Independent log-coordinate quadrature agrees with the whole cutoff
formula, including its regular missing-width term, to 1.1e−11
absolute over the tested guard ladder.

The reported finite part uses the declared radial scale 3e−4.
The identity is only for these two bare polar-delta contributions
on exact Born support. It does not include the other cuts, their
common statistical prefactor, dispersive partners, dressed lines,
or a physical damping rate.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_shared_polar_delta_reference.py
```
[Raw two-frame and radial receipt](results/finite-mass-qed-shared-polar-delta-2026-11-27.json).
