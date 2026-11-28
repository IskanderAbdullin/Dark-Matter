# Resolve the same-charge interference in the passive warm exposure

The [positive passive exposure](INDEPENDENT-WARM-BALLISTIC-HEAVY-EXPOSURE-2026-11-26.md)
bounded the same-charge reactions by Cauchy. On the same
unchanged redshifting spectra, replace that bound by the
unpolarized vacuum t/u interference. For the two equal-charge stable channels, the
two identical final heavy fermions give a 1/2 phase-space
factor. Equal populations in the two charges and the two
identical initial-pair factors sum to one relative
initial-pair normalization.

For the t/u exchange `A_t−A_u`, the mixed spin trace is
`Tr[(p₃ slash+M)(p₁ slash+m)(p₄ slash+M)(p₂ slash+m)]/4`.
Its numerator in the center of mass is
`C₀ + D² cos²θ/2`, where `D=2 p_in p_out`; the virtual product
is `A²−D² cos²θ`. Integrate that rational function analytically.
A separate 180-node log-transfer integration of the direct
four-matrix trace recovers the same-charge angular function
to 7.21×10⁻¹³ over the checked invariant masses. Near threshold
the same-charge sum is half the opposite-charge rate, not
twice it; interference need not reduce the rate by this factor
at all energies. The opposite-charge folding reproduces the
archived result at each matched quadrature to better than
2×10⁻¹³ relative.

| Anchor shape | Exact same-/opposite-charge exposure ratio | All explicitly computed channels, ΔN/N₁ | Fraction of earlier Cauchy-credit exposure |
| --- | ---: | ---: | ---: |
| q16 | .59082 | 4.55528×10⁻⁷ | .53777 |
| q32 | .59703 | 5.83594×10⁻⁷ | .54053 |
| q32 cropped to q16 | .59086 | 4.61111×10⁻⁷ | .53776 |

The q32 fixed-state crop removes about 21% of the new passive
exposure, comparable to the earlier Cauchy diagnostic.
Final Pauli factors are omitted and the formed parents are
credited with a net particle each. Only the specified passive
vacuum formation kernels are calculated exactly in angle;
this is not a full interacting evolution with inverse reactions,
inflow or repopulation. The broader [isolated energy envelope](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
continues to use a conservative all-momentum bound.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_same_charge_heavy_recycling_reference.py EOS.dat \
  "$D"/warm-ballistic-heavy-exposure-2026-11-26.json \
  "$D"/warm-q{16,32}-10p442-heavy-input-shape-2026-11-26.json
~~~
