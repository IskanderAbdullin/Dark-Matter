# Pair sensitivity to thermal replacements and extra exterior particles

## Replace the retained shape without changing number or energy

At the actual q32 7.394049 MeV endpoint and the common-chart q64
7.754805 MeV replay, independently fit one BE shape for scalar and
one FD shape for light. Use the same finite momentum domain as
the original distributions and match each represented species'
number and total energy without using its shell shape. The independent
sum agrees below 6e-12 relative. The fitted scalar/light temperatures
are 6.805/6.710 MeV on q32 and 7.132/7.052 MeV on q64. This procedure
does not impose chemical or thermal equality.

On these two examples the scalar shell-number L1 discrepancy is
1.43%, while the light discrepancy is 0.49% and 0.32%. On the
finer common collision rule, the frozen light net number rate divided
by its own N H goes from −0.34218 to −0.33990 at q32 and
from −0.33731 to −0.33530 at q64 after the replacement.
The net energy rate divided by the total retained dark density times
H changes by 1.17% and 1.08%. The changes are much larger than
the two collision quadratures' spread for this frozen observable.
The selected net pair direction is preserved in both test cases.
This is a tested moment replacement, not a proof that number and
energy alone constrain arbitrary distribution shapes.

## Conditions that allow the exterior to be nonempty

For the [frozen half-filled tagged-mode sign](INDEPENDENT-WARM-FROZEN-PAIR-BARRIER-2026-11-26.md),
any exterior **fermion** occupation in the Pauli interval [0,1]
can only decrease gain and increase forward loss; vacuum is already
the worst case for that contribution. Additional exterior **scalar**
occupation is different. Suppose its incoming support is below a
declared finite momentum p_max and its per-mode value is strictly
below a common supremum β. Positive new inverse events involve
either one exterior scalar with product at most β B_max
or both exterior scalars with product at most β², where B_max
is the maximum occupied retained scalar mode. These subsets are
disjoint. Bound their extra gain by the maximum of these two
products times the full allowed phase volume and the conservative
Yukawa amplitude envelope. Discard the extra forward stimulation
when forming this upper bound.

For p_max=0.60 GeV the sufficient β values are 8.84e-14 and
1.98e-14 on the same respective frozen q32 and q64 states.
For p_max=1 GeV the sufficient β values tighten to 4.18e-15 and
9.37e-16. These are conservative sufficient conditions for the
selected scalar-pair sign at f=1/2, not estimates of actual
exterior occupations. The stored moments do not establish these
supremum conditions, and a larger supremum does not itself prove
a reversal. No statements about prompt, elastic or a live cap follow.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
q32=results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz
q64=results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz
python -W error warm_pair_same_moment_thermal_fit_reference.py EOS.dat "$q32" "$q64"
python -W error warm_pair_exterior_scalar_allowance_reference.py "$q32" "$q64"
~~~

[Same-moment shape and rate controls](results/warm-pair-same-moment-thermal-fit-2026-11-26.json);
[exterior scalar allowances](results/warm-pair-exterior-scalar-allowance-2026-11-26.json).
