# Bounded hard tail through M/T=14: separate Born and strict extensions

The direct predecessor is the [strict x=12 result](INDEPENDENT-HARD-VECTOR-NLO-X12-2026-11-26.md).
Here the default x=12 calculation is unchanged. The extended variant
first widens the **massless Born parent** to x=14 while keeping the
strict increment at x=12 as a control. Only then is the strict
correction extended to x=14, on the complete x=13 and x=14 native
strips (the same 16 boost nodes per strip). A strict input exceeding
the parent Born interval is explicitly rejected; x=16 is not used.
The two vacuum measured bands remain disjoint from both thermal pieces.

At αs=0.30 with PCHIP/order 8, the optimistic count roots are:

| Lower hard T [MeV] | Born 12 + strict 12 [MeV] | Born 14 + strict 12 [MeV] | Born 14 + strict 14 [MeV] |
| ---: | ---: | ---: | ---: |
| 250 | 323.059 | 322.426 | 322.363 |
| 300 | 336.303 | 335.819 | 335.771 |

The Born-only extension in the middle column changes these
conditional roots by -0.633 and -0.484 MeV. At fixed expanded
parent, adding strict support changes them by only -0.063 and
-0.048 MeV, respectively. The strict PCHIP extension also changes
the slope in the previous endpoint interval; the independent linear
test isolates just the added x∈[12,14] strip. Both shifts are finite
scope comparisons, **not** estimates of a full NLO/QCD uncertainty
or proof of convergence to a physical reheating temperature.

At the tested points around each mass-window crossing, direct
physical invariant-mass quadrature verifies the isolated added
linear strip to within 4.1e-13 of the full strict increment.
Independently over the same geometry, the unsplit Born kernel checks
the repartition to within 5.6e-16, and the tested cumulative
Born order-8 versus order-12 differences are below 7.3e-12.
At the eight separate midpoints the largest relative strict-coefficient
interpolation differences are 0.027% (PCHIP) and 0.060% (linear).
None of these checks quantifies absent QCD physics or other thermal regions.
The added support above the top measured band begins at
3.671/14 = 262.214 MeV for the top x=14 edge; the lower edge
x=12 reaches it only at 305.917 MeV.

### Further Born-only range diagnostic

With no strict correction and an explicit visible clock cap of
345 MeV, compare disjoint Born edges 14, 16, 20, 24. The nominal
maximum included invariant mass is only 8.28 GeV, below the 10 GeV
vector pole; the pole is not part of this diagnostic. At the 250 MeV
lower cut, the optimistic Born roots are 324.498, 324.202, 324.079,
and 324.069 MeV. At the 300 MeV lower cut, they are 336.791, 336.616,
336.542, and 336.535 MeV. The allowed *additional Born support*
still moves the conditional benchmark more than the strict x=12–14
increment. This scope comparison is not a bound on contributions
outside x=24, on the vector resonance, or on non-Born physics.

## Reproduction

Use the original EOS, KLOE and BESIII inputs recorded in the
preceding analysis. Build the native expansion with the same
strict trace probe and 16 boost coordinates as the old strips.
The added midpoint controls are separate native runs at x=12.5
and 13.5 (u=0.25,0.7,2,5), supplied as checks, not fit points.

~~~sh
export OPENBLAS_NUM_THREADS=1
A=results/hard-vector-nlo-through-x12-strips-2026-11-26.txt
B=results/hard-vector-nlo-through-x14-strips-2026-11-26.txt
C=results/hard-vector-nlo-through-x14-checks-2026-11-26.txt
P=results/partial-hard-vector-nlo-through-x12-reheat-2026-11-26.json
Q=results/partial-hard-vector-parent14-strict12-2026-11-26.json
R=results/partial-hard-vector-nlo-through-x14-reheat-2026-11-26.json
python -W error partial_hard_vector_nlo_reheat_reference.py EOS.dat KLOExSec.data BESIII-2112.11728.txt "$A" --born-x-max 14 >"$Q"
python -W error partial_hard_vector_nlo_reheat_reference.py EOS.dat KLOExSec.data BESIII-2112.11728.txt "$B" --born-x-max 14 --intermediate-check "$C" >"$R"
python -W error hard_vector_thermal_extension_audit_reference.py EOS.dat "$A" "$B" "$P" "$Q" "$R"
python -W error hard_vector_nlo_overlap_audit_reference.py EOS.dat "$B" --born-x-max 14
python -W error hard_vector_born_subpole_scope_reference.py EOS.dat KLOExSec.data BESIII-2112.11728.txt
~~~
