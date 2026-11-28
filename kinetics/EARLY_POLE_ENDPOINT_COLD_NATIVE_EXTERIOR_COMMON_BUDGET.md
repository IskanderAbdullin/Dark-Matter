# A common exterior does not by itself erase the retained cold moment ordering

Use the two accepted sign-targeted windows and normalize each whole
response separately. Suppose a common exterior component has exactly
the same normalized mass h and the same conditional distribution in
both charts. That common component cancels, regardless of its location
or within-component shape, provided the moment exists. It is the
unmatched exterior that can reverse the retained ordering.

Let a be the smallest compatible retained fine moment, b the largest
retained third moment, and H the upper exterior probe. The unmatched
fractions eta_fine and eta_third can keep the original sign for all
compatible shapes exactly when

    a * eta_fine + (H - b) * eta_third < (1 - h) * (a - b).

The adverse independent allocations put the unmatched fine response at
zero and unmatched third response at the upper cap. The same reasoning
applies to RMS after working with the second moment. The thresholds
are sharp at the closure of the support.

For forward-number mean q with formal cap .0015 and no common exterior
mass, the limiting unmatched fraction if only the fine chart has one
is 0.00368774. If only the third has one it is 9.41128e−6. If both
have equally much independently shaped unmatched response the limit is
9.38732e−6, as in the earlier independent result. Equal *fully shared*
exterior contributions cannot reverse the order by themselves.
No shared or unmatched fraction has been measured here.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_cold_native_exterior_common_budget.py --root results \
 --sign results/early-pole-cold-native-moment-sign-refinement-2026-11-28.json
```
