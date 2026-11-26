# Reconstruct the discrete cold-cell budget on the corrected q64 trajectory

For four intervals of the
[corrected q64 trajectory](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md),
rerun the actual accepted split step with the identical parent,
clock, archived exterior gravitation, and q32 Newton initial guess.
All four replay snapshots have the exact stored hashes. Track the
smallest light cell through both redshifts, both conservative pair
kicks and the implicit prompt-plus-below-pole elastic update.

Use the two actual **midpoint** pair rates for the kicks, not the
frozen endpoint pair rate. The centered stiff residual is
backward Euler: the local prompt Jacobian and the separate
elastic evaluator supply the two occupation source terms at
the accepted implicit state. Their sum matches the combined
evaluator. The recombined cell changes close to below 1e-16 in
occupation; each redshift changes the comoving cell average by
less than 1e-17. A collisionless redshift is not a competing
loss in this comoving observable.

| q64 interval end | Pair per ln(a) | Prompt per ln(a) | Below-pole elastic per ln(a) | Total |
| ---: | ---: | ---: | ---: | ---: |
| 7.941700 MeV | −1.10246e-5 | −5.37171e-8 | −3.81089e-5 | −4.91873e-5 |
| 7.754805 MeV | −1.07556e-5 | −2.60446e-8 | −3.92455e-5 | −5.00272e-5 |
| 7.572290 MeV | −1.04487e-5 | −1.21217e-8 | −4.00505e-5 | −5.05113e-5 |
| 7.394049 MeV | −1.01106e-5 | −5.42827e-9 | −4.06667e-5 | −5.07827e-5 |

Elastic redistribution accounts for about 77–80% of the tested
cold-cell decrease, with the pair kicks contributing most of
the remainder. The prompt contribution is negative and small.
This is attribution in a particular discrete split scheme,
not a claim about independent runs with each interaction off.

At the final 7.394049 MeV q32-48 and q64-24 snapshots, perform
a separate frozen first-six-and-nine-cell angular/radial sensitivity
check. The nine-cell tagged rate is negative on all nine cold cells
at each of three quadrature orders. The largest
relative rate difference between successive rules is 3.58e-9.
This last test leaves the native complement and elastic operator
unchanged. It does not bound unrepresented phase, interpolation,
or complete physical angular error.

Use the local `ladder/` workspaces reproduced in the linked note.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for j in 6 12 18 24; do
  python -W error warm_q64_corrected_actual_cold_budget_reference.py EOS.dat \
    "ladder/wide/q64/step-$((j-1)).json" "ladder/wide/q64/step-$j.json" \
    "ladder/q32-24/step-$j.json"
done
python -W error warm_cold_pair_multishell_audit_reference.py EOS.dat \
  results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
  results/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz --finer
~~~

[Four full split budgets](results/warm-q64-corrected-actual-cold-budget-2026-11-26.json);
[frozen terminal tagged quadratures](results/warm-corrected-terminal-cold-angular-refinement-2026-11-26.json.gz).
