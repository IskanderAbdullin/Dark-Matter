# Translate the remaining 50 MeV source budget into spectral units

Use the 50 MeV [charged two-pion control](INDEPENDENT-WARM-PION-STATISTICS-AND-THRESHOLD-2026-11-26.md)
with the wider conditional form-factor caps as the baseline.
Consider an additional positive source, described only by
an **effective** multiplier `R_X(T,s)` of the unblocked
unit electromagnetic spectral kernel. This definition
puts any applicable population factors in `R_X`; it does
not identify the multiplier with a measured vacuum `R`.
Use charge-neutral portal pair injection with up to one
direct count per m₂ of energy, then credit all the extra
comoving energy at the 50 MeV start. The extra input needed
for this sufficient comparison to reach .060 is about
1.082×10⁻⁹ GeV⁴ in start-normalized units.

If `R_X(T,s) ≤ C_j` over each disjoint band and time, the total
extra source energy is at most `Σ_j C_j E_j`. Bound each
`E_j` with positive invariant and visible-temperature cells,
the sub-pole propagator below 2 GeV, and the analytic width floor
and Bessel upper tails above it. The interval below .28 GeV is
an illustrative scan for a low-mass **unmeasured** input,
not an additional two-charged-pion channel. There is no
double assignment of the known source: `R_X` is additional.

| Extra invariant √s band, GeV | Upper energy for unit C, GeV⁴ | Cap C touching .060 if this is the only extra band |
| --- | ---: | ---: |
| .130–.280 | 5.804×10⁻¹⁰ | 1.86 |
| .280–.316 | 1.833×10⁻¹⁰ | 5.90 |
| .316–.420 | 3.789×10⁻¹⁰ | 2.86 |
| .420–.600 | 2.216×10⁻¹⁰ | 4.88 |
| .600–.975 | 3.473×10⁻¹¹ | 31.2 |
| .975–2.0 | 2.098×10⁻¹³ | 5.16×10³ |
| 2.0–∞ | 3.778×10⁻¹⁴ | 2.86×10⁴ |

Low-invariant coverage is more important per unit spectral
strength. The individual caps are **not simultaneous**;
when multiple bands are added, compare their sum with the
single comoving budget. None of these caps is asserted for
a real hadronic medium. The relatively loose width-floor
bound above 2 GeV is not an estimate of the pole rate.
Additional heat, inflow of externally made light scalars
and sources below .130 GeV need separate allowances.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_unmeasured_spectral_budget_reference.py EOS.dat \
  "$D"/warm-pion-bose-and-threshold-2026-11-26.json
~~~
