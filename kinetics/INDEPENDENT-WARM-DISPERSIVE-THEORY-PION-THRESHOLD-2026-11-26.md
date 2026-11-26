# Bound the sub-0.30 GeV charged-pion slice with a separate theory input

The [BABAR transfer](INDEPENDENT-WARM-BABAR-DISPERSIVE-GAP-2026-11-26.md)
does not constrain charged pion pairs between the physical
threshold and .30 GeV. Use a distinct primary theory result:
[Ananthanarayan, Caprini and Das, Eq. (23)](https://arxiv.org/pdf/1810.09265)
report the positive dispersive moment up to .63 GeV as
`(132.91 ± .76)×10⁻¹⁰` with the electron-positron inputs;
it uses a form-factor determination based on analyticity
and unitarity. The paper includes the one-photon final-state
radiation correction and excludes vacuum polarization.
This is a **conditional theory prior**, not a direct measurement
below .30, and its quoted uncertainty is used only to choose
a numerical input allowance.

Reallocate the entire moment to the interval `2mπ≤√s≤.30`
using the integrated thermal upper/lower dispersive kernel
ratio. This strongly overestimates a model constrained by
actual threshold shape, but avoids interpolation assumptions.
Use independent intervals: theory only below .30, the BABAR
allowance only .30–√.10 GeV, and the earlier Bose-weighted
KLOE effective band above √.10. Do not combine these into a
joint dispersion estimate or assume statistical agreement.

The refined threshold ratio is at most .01691161 GeV⁴
(coarser .01747145). With the theory central value plus
twice the quoted error, the extra threshold-slice energy
upper is 2.27343×10⁻¹⁰ GeV⁴; the separate BABAR upper is
1.21881×10⁻¹⁰. These two allowances together spend 31.55%
of the previously remaining source budget. A conservative
convex interpolation then gives Ωφh²≤.03826 if nothing
else is added. The remaining input budget is
7.5783×10⁻¹⁰ GeV⁴ under the same positive early-deposit
accounting. The wider pion elastic and neutral pion-photon
tests remain separate from this comparison.

No empirical low-bin pion spectrum, full finite-temperature
source or relic prediction follows from the result. It retains
the effective vacuum two-pion assumption and the earlier
conditional normalization, including Bose population factors.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_dispersive_theory_pion_threshold_reference.py EOS.dat \
  "$D"/warm-babar-dispersive-gap-2026-11-26.json
~~~
