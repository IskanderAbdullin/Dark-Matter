# How large would a short Bose burst have to be?

The [fermion floor with vacuum lepton loss](INDEPENDENT-WARM-FERMION-SURVIVAL-WITH-LEPTONS-2026-11-26.md)
used one outgoing-occupation cap for all time. The
negative quadratic internal loss comparison also works
for a time-varying cap `F(b)`. In its analytic Riccati
solution replace the internal factor `I(1+F)^2` by
`I_B=∫ b⁻³(1+F(b))²/H d ln b`. The vacuum lepton
quadratic contribution uses its old integral I without
Bose factors; the separate high-invariant energy loss
credit is unchanged. For the q32 snapshot the target
.120 cannot be reached by this network unless this
weighted Bose integral is at least 34.75 times the
previous nominal **upper** expansion integral.

Test a declared outside cap `F_out` and a larger
inside cap `F_in` over one burst window. With J the
upper expansion integral on that window, the
conservative input to the comparison is
`I_upper(1+F_out)^2 + J[(1+F_in)^2−(1+F_out)^2]`.
The upper integrals are bounded at the visible-EOS
knots, so no arbitrary photon-clock interpolation
is required. The ratios J/I_upper below compare
two upper integrals; they are not measured time
fractions of an interacting trajectory.

| Relative expansion window after 10.442 MeV | J/I_upper | Necessary allowed F_in if F_out≤1 | If F_out≤.1 |
| --- | ---: | ---: | ---: |
| 1–1.12 | .10704 | 16.07 | 16.74 |
| 1–1.5 | .33297 | 8.82 | 9.10 |
| 1–3 | .66589 | 6.08 | 6.18 |
| 1.5–3 | .33292 | 8.82 | 9.10 |
| 3–10 | .23294 | 10.66 | 11.05 |

A cap lower than the listed one still rules out the
target total in this limited comparison. These are
alternatives, not simultaneous demands for overlapping
windows. They say where a putative numerical solution
would have to find exceptionally large outgoing-state
occupancy or introduce another loss mechanism. They do
not infer a physical Bose burst, establish either cap,
or apply to additional hadronic/visible loss channels.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_bose_weighted_depletion_budget_reference.py EOS.dat \
 "$D"/warm-isolated-fermion-survival-bound-2026-11-26.json \
 "$D"/warm-fermion-survival-with-leptons-2026-11-26.json \
 "$D"/warm-isolated-recycling-energy-envelope-2026-11-26.json
~~~
