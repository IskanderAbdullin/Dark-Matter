# Match real and virtual photons in an isotropic low-temperature rest bath

Use a photon bath at rest with the parent, equal fermion charges,
and T much less than both charged masses and the available decay
energy. The archived spectra supply two mass ratios only; charged
occupations and final-state occupation factors are deliberately absent.
This does not use the archive's physical temperature or its
relativistically moving parent distribution.

The benchmark is Section IV of [Czarnecki et al.](https://arxiv.org/abs/1110.2171),
which treats a small daughter mass, a massless neutral scalar and a
chiral current. Their one-loop photon pole shift and spin sum
(Eqs. 24, 30, 32–33) are also applicable to the external fermions
in the scalar-current control. Write h = [g²/(4 pi²)] (pi²/6) T².
Then the squared pole masses of both charged legs shift by +4h
while the spin numerator contains the **vacuum** mass and the
combination slash(P−K). The shifted two-body phase space must be
varied together with these spin sums. The finite projection of an
isolated fixed-spin derivative is not another term to add on top of
this pole-expansion convention.

At each tested photon direction, the leading signed real pair,
the null-photon vertex interference and the singular pole residues
cancel within 5.4e−16 relative. The null-photon vertex has only this
leading term after pairing its two signs. For the real remainder,
evaluate the two signed recoils with a gamma-matrix current and
take the even second derivative. For the scalar current an
independent complex contour of the reduced invariant amplitude
gives the same derivative; both currents also have two independent
angular-node and finite-difference sequences. The remaining real
and external terms are individually sensitive to the masses.

On the fine mass fixture the finite real coefficient is +18.78012034;
the external phase-space and spin-sum response gives −20.78012033.
The sum is −2 within 6.3e−9. The third mass fixture agrees.
All 18 matched controls, including sizable daughter and scalar
masses and both chiral and scalar currents, agree with −2 within
3.2e−8. This extends the numerical check beyond the published
small-mass fixture; it is not an analytic universality proof.
The correction is therefore consistent in these controls with
−g²(T/M)²/12 relative to Born. Equivalently the coefficient −2
multiplies [g²/(4 pi²)] (pi²/6) (T/M)².

Without replacing the real integrand by its leading coefficient,
the direct finite-bath ratios at T/M=.001, .002, .004 are approximately
−2.0000174, −2.0000695 and −2.0002782. The deviation from the leading
term is quadratic in T/M after dividing by (T/M)², as expected
for the next even recoil term.

No charged thermal cuts, plasma screening, finite widths, final-state
Bose/Pauli weights, boosted bath, or self-consistent cosmological
evolution is included. In particular this reference does not justify
extrapolation to the high-temperature production epoch.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_matched_low_temperature_photon_bath_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-matched-low-temperature-photon-bath-2026-11-27.json).
