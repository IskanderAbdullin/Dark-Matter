# A leading-soft subtraction for the charged Yukawa decay

The [finite-cut real radiative test](INDEPENDENT-EARLY-POLE-YUKAWA-REAL-RADIATIVE-WITNESS-2026-11-27.md)
does not constrain radiation below its photon cut. Before extending
it, a leading-soft real/virtual check is needed. The soft branch in
this reference is factored from a fixed Born event: equal charges,
a neutral emitted scalar, and massive on-shell fermions.
The thermal real/virtual pattern is the leading eikonal one of
[Weldon, 1993](https://arxiv.org/abs/hep-ph/9311245), not a finite
one-loop virtual calculation for this decay.

For photon direction n=(1,n̂), the current is
`j = P/(P·n) − b/(b·n)`, and n·j=0. With relative rapidity
cosh(χ)=P·b/(M m), the angular average of −j² is
`B=2(χ coth χ−1)`. The photon angular integral and phase-space
normalization multiply dω/ω by g²B/(4π²) with the amplitude
normalization of this check. In the parent rest frame the weight
is β²(1−z²)/(1−βz)². The independently evaluated plasma-frame
integral resolves the narrow forward structure after analytic
azimuthal integration. Both archived scalar pole kinematics at
Q=4.1e−6 give B≈3.6669, matching the direct rest and plasma angular
integrals to 6e−13 absolute. Comparing the coherent exact real
amplitude at six rest-frame angles shows linear approach to
the soft limit, with relative error at most 3.0e−9 for ω/M=1e−9.

For a measurement F under signed soft energy shifts, combine
three leading terms **before** integrating:
`(1+n_B) F(+ω)+n_B F(−ω)−(1+2n_B) F(0)`.
The thermal part is a symmetric second difference, and the
vacuum part is a first difference. For a smooth F the combination
is O(ω) even though the real and counterterm pieces separately
grow as 1/infrared cutoff after integration. An exponential test
weight agrees with an independently combined 80-digit integral
to 7e−18, and its integral approaches the zero-cut result
linearly. An intentional `F=1+c|ω|` cusp instead leaves a
logarithm with coefficient 2cT, reproduced by the cutoff ladder.
This manufactured cusp is a warning about regularity of the
inclusive event measure; it does not show that an actual
integrated archived shell observable has this infrared behavior.

This calculation fixes the universal leading soft factor and tests
the subtraction structure only. Finite virtual terms, angle-dependent
recoil of the real operators, resummed in-medium photon propagators,
and matching to the existing finite-cut real channels remain open.
No physical radiative rate or upper bound follows.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_soft_subtraction_reference.py
```
[Raw receipt](results/finite-mass-qed-soft-subtraction-2026-11-27.json).
