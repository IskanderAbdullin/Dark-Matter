# Neutral-pion photon input below 0.60 GeV with a conditional form-factor cap

One low-invariant possibility left open by the
[spectral budget](INDEPENDENT-WARM-UNMEASURED-SPECTRAL-BUDGET-2026-11-26.md)
is vacuum `π⁰γ` annihilation. Adopt the form-factor convention
of [Hoferichter et al., Eq. (34)](https://doi.org/10.1140/epjc/s10052-014-3180-0)
and define `Fhat=F/e²`, so the spectral weight in units of
the point-muon cross section is
`R=(π α/2) s (1−mπ⁰²/s)^3 |Fhat|²`.
No form factor is extrapolated from the [SND data](https://arxiv.org/abs/1601.08061),
whose reported coverage starts at 0.60 GeV.

For a massless initial photon, a uniform massive-particle
Bose multiplier would be invalid. At fixed `s` and total
lab energy `E_sum`, even the most deboosted photon satisfies
`Eγ ≥ (s−mπ⁰²)/(4 E_sum)`. Thus its population is bounded
by its MB weight times `1+4 T E_sum/(s−mπ⁰²)`.
Use the separate massive pion factor for a zero-chemical-
potential bath. The rate energy moment now uses the first
and second thermal moments. The cubic spectral suppression
cancels the apparent threshold denominator with two powers
to spare. Positive cells bound both terms, including a
massless analytic time tail at b≥500.

| Covered √s from mπ⁰ to | Upper comoving energy per unit (GeV⁻¹)² of squared Fhat | Largest uniform Fhat magnitude that fits the remaining upper budget |
| --- | ---: | ---: |
| .28 GeV | 2.2510×10⁻¹³ GeV⁴ | 69.34 GeV⁻¹ |
| .42 GeV | 1.0562×10⁻¹² GeV⁴ | 32.01 GeV⁻¹ |
| .60 GeV | 1.8068×10⁻¹² GeV⁴ | 24.47 GeV⁻¹ |

These windows are nested, not additive. The ceilings are
sufficient conditions for this one added channel under the
chosen vacuum kernel and clock. They are **not** estimated
form-factor magnitudes or a bound on other sources, and no
thermal spectral function or plasma photon model is inferred.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_neutral_pion_photon_budget_reference.py EOS.dat \
  "$D"/warm-unmeasured-spectral-budget-2026-11-26.json
~~~
