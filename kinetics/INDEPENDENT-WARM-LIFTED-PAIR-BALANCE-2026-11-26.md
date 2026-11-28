# Detailed-balance and late first-shell controls for the cold lift

Test the [nine-cell reverse partition](INDEPENDENT-WARM-COLD-PAIR-MULTISHELL-2026-11-26.md)
on a synthetic FD/BE bath at a common 15 MeV temperature and common
−20 MeV chemical potential. Over the interior of the finite momentum
domain the gain and loss integrands must cancel; occupations above
the represented edge are set to zero and physical final-state vacuum
is still allowed. The check is specifically on the nine very narrow
fermion cells, where the distance from the cutoff makes that thermal
tail truncation immaterial at the observed precision.

Normalize each cell's net rate by its positive forward loss. The
maximum absolute relative defect is 2.21, 0.253 and 0.371 for the
three original quadratures. Merely increasing the original angular
order on this grid does not reliably remove the failed cancellation.
The corresponding nine-cell defects are 6.86e-7, 7.42e-10 and
3.91e-12, keeping the forward rule unchanged in each comparison.
The larger residuals follow the order of the still-discrete forward
integral. An independent 23 MeV bath, with chemical potential −28 MeV and
a 0.90 GeV outer edge, gives 0.414 native versus 1.35e-12 lifted at the
finer order. This is a local detailed-balance test, not a guarantee of
cancellation in the unmodified inverse sector.

Independently freeze the actual corrected q32 7.394049 MeV
[endpoint](INDEPENDENT-WARM-MULTISHELL-CORRECTED-TRAJECTORY-2026-11-26.md).
With the finer below-pole elastic control, its signed first-shell
rates are zero for the selected prompt, −4.0987e-5 for elastic and
−1.0082e-5 for the conservative scalar pair contribution, per ln(a).
The total −5.1069e-5 is close to the actual last-step decline
near −5.08e-5. This confirms that the decline persists after removal
of the native pair-rate aliasing and that the elastic term dominates
this local selected budget. The result does not test unselected
interactions or a whole-spectrum elastic convergence ladder.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_multishell_pair_balance_reference.py
python -W error warm_multishell_trajectory_local_budget_reference.py EOS.dat \
 results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz
~~~

[Per-shell balance residuals](results/warm-multishell-thermal-balance-2026-11-26.json)
and [frozen actual endpoint rates](results/warm-multishell-local-budget-q32-endpoint-2026-11-26.json).
