# A frozen stiff collision-step integrator comparison

The [q16 warm refinement](INDEPENDENT-WARM-EARLY-TIME-REFINEMENT-2026-11-26.md)
found a substantial first-order component in the evolving operator
split. Here isolate **one fixed 12 MeV prompt-plus-elastic collision
problem**, with the same two physical moments, support, collision
panels and centered fluxes. There is no FRW remapping, pair kick,
new source, or evolution of collision geometry in this comparison.

The original solver uses backward Euler. To test implicit midpoint,
solve the half-step moment equation using the same centered collision
flux and then reflect its moments about the initial value. For the
two-stage L-stable alternative (SDIRK2), let
`g = 1 - 1/sqrt(2)`. The two implicit equations are

~~~text
U1 - U0 = g h F(U1)
U2 - [U0 + (1-g) h F(U1)] = g h F(U2).
~~~

The bracket is a moment **right-hand side**, not a state whose
occupation must be reconstructed. Both actual stages still use
the original positive entropy variables. An explicitly supplied
offset changes only the moment target of the Newton equations;
it leaves the old default and its flux unchanged. A legacy warm
bridge with no offset is reproduced byte for byte. The optional
two-stage method is not unconditionally positivity-preserving.

## Tested behavior

With the same total frozen normalized collision interval 0.002,
compare 2, 4 and 8 substeps for each method. All successful endpoints
conserve species number and total energy and increase total entropy
in this run. For the largest of the four species/shell-moment
weighted differences across 4→8 steps:

| Method | Max. relative shell-moment change |
| --- | ---: |
| Backward Euler | 2.16×10⁻⁶ |
| Implicit midpoint | 1.25×10⁻³ |
| L-stable two-stage method | 4.74×10⁻⁷ |

The midpoint change is actually worse than its own 2→4 difference.
This is consistent with its undamped limit for extremely stiff
linear modes, not a measurement that the formal order is negative.
For a linear mode the two-stage method has stability function
`(1 + (1-2g)z)/(1-gz)^2`, which tends to zero for large negative z;
midpoint tends to −1. The varied empirical changes are not used to
claim a convergence order in this stiff pre-asymptotic test.

This is a numerical control on one frozen operator and is
**not the default expanding warm method**. A larger step or another
support may challenge realizability, conservation tolerances, or
positivity. The result does not resolve the full warm time error
or establish faster whole-run performance.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error warm_stiff_microstep_refinement_reference.py \
  EOS.dat results/50MeV-from-empty-pair-quantum-shape-2026-11-26.json \
  --step .002 --steps 2 4 8
~~~
