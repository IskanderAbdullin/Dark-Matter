# Near-rest Bose enhancement can change the omitted-source scaling

The [rest-strip bound](INDEPENDENT-TWO-MASS-REST-STRIP-EXACT-AREA-2026-11-27.md)
controls unweighted geometry. To see what changes once occupations
are added, use thermal Bose and Fermi distributions for a positive
on-shell mass pair, with fully signed forward and inverse
Bose/Pauli factors. Keep the kinetic offsets, rest contributions,
and (possibly different) inverse temperatures separate.

For the decay and inverse reaction, write A=log(F/R). Avoid subtracting
the logarithms of individually large flows or the daughter energies.
The thermal formula is

    A = A_rest + (1/T_s−1/T_P) K_s + (1/T_f−1/T_P) K_f.

Here A_rest includes the separately supplied rest-energy and chemical
terms. With exact common-temperature chemical balance A_rest=0.
Reconstruct F−R with expm1 at fixed affinity and keep the original Bose
exponent strictly positive on every sampled point. Then the output
has the sign of A and event entropy production (F−R)A is nonnegative.
At exact balance the event source is zero. For a near-threshold
event with |A|=10⁻²⁸, conventional float64 forward-minus-reverse
subtraction gives zero while the stable signed source is ±3.66252e−8
in the fixture units. Independent 100-digit direct quadrature-point
evaluation checks these small signed sources; they are not rates.

## Three different rest-strip limits

With p the boson momentum, its energy offset K_s and a strictly
positive Bose exponent, |F−R|≤(1+f_B) min(1,|A|). The unweighted
dE_s dE_f area scales cubically in the upper momentum cutoff p_c.

- **Gapped occupation.** An exponent floor δ_B>0 bounds 1+f_B by
  1/(1−exp(−δ_B)), preserving the cubic bound.
- **Saturated thermal occupation and nonzero rest affinity.**
  At p>0, 1+f_B≤1+T_s/K_s. The exact extra strip weight is
  ∫dA/K_s=(2D/m²)[p_c+m asinh(p_c/m)], so the bound is linear
  rather than cubic. The p=0 point is never evaluated.
- **Saturated occupation with tuned affinity.** If only the scalar
  temperature differs and A_rest=0, the event affinity is proportional
  to K_s. This cancels the infrared factor and the cubic law returns.

The manufactured strip integrals obey these separate bounds.
Two quadrature orders recover exponents 1, 3 and 3, respectively,
over the smallest tested cutoffs. Therefore the unweighted area
alone cannot certify a cubic error for a weighted source. The
tests do not assume a physical occupation, thermal history, or
effective rate for this project.

Replay from the kinetics directory with NumPy, SciPy and mpmath:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_quantum_kinetic_affinity_reference.py
```
[Full result](results/two-mass-quantum-kinetic-affinity-2026-11-27.json).
