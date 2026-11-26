# Extend the controlled hard strips through M/T = 16

This follows the [separate x=14 Born and strict extensions](INDEPENDENT-HARD-VECTOR-NLO-X14-2026-11-26.md).
The formerly missing x=15 strip was explicitly calculated on all
16 boost nodes before using the available x=16 strip. Separately
first extend the disjoint massless Born parent to x=16, then apply
the x=14–16 strict increment **within that same parent**. No new
interpolation is performed across the remaining x=16 to x=20 gap.
All strict strips still exclude both previously measured vacuum
bands. The default x=12 case and its input/output formats are
not changed.

At illustrative αs=0.30 (PCHIP/order 8), the conditional
all-converted-count roots are:

| Lower visible T [MeV] | Born 14 + strict 14 [MeV] | Born 16 + strict 14 [MeV] | Born 16 + strict 16 [MeV] |
| ---: | ---: | ---: | ---: |
| 250 | 322.363 | 322.084 | 322.057 |
| 300 | 335.771 | 335.605 | 335.588 |

With the broader parent fixed, the strict extension moves the
root by −0.0273 and −0.0163 MeV. These values include PCHIP's
small slope change on the old endpoint panel. An independent
physical-mass calculation isolates only the new strips in
the linear case and agrees within 1.8×10⁻¹³ relative to the
whole local strict rate near the marked edges. The comparison
between the linear/order 4 and PCHIP/order 8 control at αs=.3
is ≤0.00046 MeV at the two roots. These are numerical tests
at fixed coupling and clock, not a QCD error estimate. The eight
held-out midpoint coefficients differ by at most 0.019% (PCHIP)
or 0.035% (linear), not a continuum guarantee.

## A separate, nonoverlapping pole proxy

The previous Born-only thermal support at x≤24 and T≤345 MeV
never reaches 8.5 GeV. Without adding anything to that result,
integrate a distinct 8.5–11.5 GeV mass window under the same
toy *three massless flavors* Born trace and 10 GeV vector
propagator. Integrating the stiff mass profile independently
in log(M) and in arctan((M²−V²)/(VΓ)) gives agreement better
than 1.9×10⁻¹⁰ at the tested settings, and order 8 versus 12
agrees to 4.8×10⁻⁷. Between 95.2% and 95.5% of this
**band's** count is in 9.95–10.05 GeV. The band alone contributes
optimistic Ωh² of 6.90×10⁻⁶ and 2.00×10⁻⁵ for the respective
250→322.36 and 300→335.77 MeV clocks. That is about
0.0058% and 0.0168% of the reference .1188 count target.

This is not measured hadronic input, a prediction of a 10 GeV
resonance contribution, or a bound on all missing regions or
heavy flavors. In particular, no gap from the earlier x cap to
8.5 GeV (nor masses above 11.5 GeV) is filled by this proxy.
It is a separate numerical sensitivity to the model weight inside
this one band, not a correction applied to either root above.

## Reproduction

Use the same primary radiation, KLOE, and BESIII files as before.
Run the native adapter on the complete x=15 strip and the separate
x=14.5/15.5 probes (four boosts each) at the same published
tolerance; see the [adapter instructions](INDEPENDENT-HARD-VECTOR-NLO-STRIPS-2026-11-26.md).
Do not use the midpoint checks as fit nodes.

~~~sh
export OPENBLAS_NUM_THREADS=1
A=results/hard-vector-nlo-through-x14-strips-2026-11-26.txt
B=results/hard-vector-nlo-through-x16-strips-2026-11-26.txt
C=results/hard-vector-nlo-through-x16-checks-2026-11-26.txt
P=results/partial-hard-vector-nlo-through-x14-reheat-2026-11-26.json
Q=results/partial-hard-vector-parent16-strict14-2026-11-26.json
R=results/partial-hard-vector-nlo-through-x16-reheat-2026-11-26.json
python -W error partial_hard_vector_nlo_reheat_reference.py EOS.dat KLOExSec.data BESIII-2112.11728.txt "$A" --born-x-max 16 >"$Q"
python -W error partial_hard_vector_nlo_reheat_reference.py EOS.dat KLOExSec.data BESIII-2112.11728.txt "$B" --born-x-max 16 --intermediate-check "$C" >"$R"
python -W error hard_vector_thermal_extension_audit_reference.py EOS.dat "$A" "$B" "$P" "$Q" "$R" --old-max 14 --new-max 16
python -W error hard_vector_nlo_overlap_audit_reference.py EOS.dat "$B" --born-x-max 16
python -W error hard_vector_born_isolated_pole_scope_reference.py EOS.dat
~~~
