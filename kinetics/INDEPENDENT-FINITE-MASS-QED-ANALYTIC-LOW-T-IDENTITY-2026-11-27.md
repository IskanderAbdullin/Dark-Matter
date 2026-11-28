# Analytic low-temperature photon-bath cancellation at finite masses

The [numerical matched benchmark](INDEPENDENT-FINITE-MASS-QED-MATCHED-LOW-TEMPERATURE-PHOTON-BATH-2026-11-27.md)
recovers the −2 coefficient in a finite set of controls. Here the
remaining leading term can be simplified symbolically for both
the scalar and chiral currents. The setup is the equal-charge
photon-only rest bath of [Czarnecki et al.](https://arxiv.org/abs/1110.2171),
extended to finite charged and neutral daughter masses in the
stated two-body channel.

Scale the parent mass to one, use charged mass m and neutral
mass s with m>0 and m+s<1, and write r for the signed real-photon
energy. Set S=1−2r, u=(1+m²−s²)/2 and p=√(u²−m²).
In the recoiling pair rest frame define
w=(E_b−p_b cosθ)/√S. Its limits satisfy

    S w² − (S + m² − s²) w + m² = 0.

The phase-space angular factor cancels against dw/dcosθ. For
the scalar current set a=m. The spin-summed chiral rate has
the same normalized expression with a=0: the parity-even
rate is the parity-even sum of the scalar and pseudoscalar spin
traces before the common Born normalization, and the spin numerator
of the pseudoscalar is obtained from
the scalar expression by reversing the charged daughter's mass
in the trace, not in the kinematics. The spin-summed scalar–pseudoscalar
interference vanishes. A direct chiral gamma-matrix replay
independently tests this reduction at signed photon energies.

The integrand after multiplying the squared amplitude by r² is
a quadratic in r at fixed w:

    C0 = −4(u+a)(1 + m²/w² − 2u/w)
    C1 =  8(u+a)(1 − 1/w)
    C2 =  4(w − 2 + 1/w).

Divide its w integral by 8p(u+a). Differentiate this integral
twice, including both moving endpoints. At r=0, w_±=u±p
and w_+ w_−=m². The real even recoil coefficient is

    G''(0) = 2 − [u(4−2u)−2]/p² − [4−2u−log(w_+/w_−)/p]/(u+a).

From the thermal pole expansion both mass squares shift by +4h,
where h=(g²/(4 pi²))(pi²/6)(T/M)² in parent-mass units. The external
phase factor and spin-sum derivative is

    C_external = [u(4−2u)−2]/p² − 4
                 + [4−2u−log(w_+/w_−)/p]/(u+a).

The mass-dependent rational and logarithmic terms cancel, leaving
G''(0)+C_external=−2 exactly. The leading infrared terms have already
cancelled between the signed real pair, vertex interference and the
pole residues in the matched reference. The symbolic test verifies
the antiderivative and the final zero residual, and separate
gamma-matrix evaluations verify the signed real formula across
emission and absorption, not just its small-r coefficient.

Thus the relative leading correction in this model is
−g²(T/M)²/12. The identity is an asymptotic result for a nondegenerate
open two-body channel; near threshold the bath must also be cold
relative to the available energy. No charged thermal occupations,
resummed photon, finite widths, moving-parent bath or hot cosmological
collision operator follows from it.

Replay with SymPy, NumPy and SciPy from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_analytic_low_t_identity_reference.py
```
[Raw symbolic and numerical receipt](results/finite-mass-qed-analytic-low-t-identity-2026-11-27.json).
