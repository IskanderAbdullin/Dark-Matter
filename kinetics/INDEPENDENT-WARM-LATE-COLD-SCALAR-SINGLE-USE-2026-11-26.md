# Limit cold scalar return from one use of the tracked exterior at 7.394 MeV

The [soft recoil gate](INDEPENDENT-WARM-LATE-COLD-SCALAR-PROMPT-GATE-2026-11-26.md)
is closed for the represented incoming spectra at p/T=.1.
This leaves open what the collisionless exterior receipts
could supply. With both inputs drawn only from the frozen
represented spectrum or its tracked exterior, each parent
above the gate must use at least one exterior particle with
energy ≥ half the gate. If exactly one partner is exterior,
a stronger threshold is the gate minus the maximum
represented energy of its partner species. For each positive
receipt, use its extrapolated birth threshold, birth floor,
mean and count to bound the eligible population; the
conditional tail-count inequality is the one used in the
[high-band receipt audit](INDEPENDENT-WARM-HEAVY-ELIGIBLE-EXTERIOR-RECEIPT-2026-11-26.md).
Each tracked exterior particle can be consumed at most
once for this diagnostic.

The purely geometric bound counts one low-momentum scalar
from every eligible heavy decay, independent of angular
statistics. It is below 7.92e−10 (q32) and 2.38e−10 (q64)
of the existing represented p/T<.1 scalar count, even
allowing two exterior inputs. No within-receipt distribution
or cutoff above the recorded floor is assumed.

## Further conditional angular allowance

For an unpolarized vacuum decay write the parent boost
rapidity as y and the scalar rest-frame rapidity as y*.
The uniform decay energy interval has width 2 p* sinh y
and minimum energy m_s cosh(y−y*). The fraction with
E_slab<E_cut is the positive difference of the cut and
this minimum divided by that width. It is positive only
within |y−y*|<asinh(p_cut/m_s), and its derivative has the
sign of E_s,rest−E_cut cosh y. The maximum is thus at
cosh y=E_s,rest/E_cut; a separate grid checks the maximum.

At p/T=.1 this best vacuum fraction is below 5.94e−5.
Multiplying the eligible one-use count by it gives
4.70e−14 and 1.42e−14 of the initial cold q32 and q64
counts, respectively. With fixed archived cold Bose
occupation and the **additional assumption** f_light≤.5
throughout the open decay phase, the fractions are bounded
by 9.40e−14 and 2.83e−14. The latter multiplies the vacuum
angular fraction by (1+max f_cold)/(1−.5), since the
total open decay weight is at least (1−.5) times its
vacuum value. The global light cap is not inferred from
the collisionless receipts.

Neither allowance controls repeated use, preexisting
parents, dynamically increased momenta or occupancy,
untracked arrivals, other collision processes, or entry
into the cold band by redshift. This is not an integrated
physical bound on the final cold scalar population.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_scalar_single_use_receipt_reference.py \
 results/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Receipt and angular allowances at three momentum cuts](results/warm-late-cold-scalar-single-use-receipts-2026-11-26.json).
