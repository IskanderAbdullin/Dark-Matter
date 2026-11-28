# Why the minimum-TV overlap survives direct subpanel splitting

The directly integrated eight-parent split made the W1 bounds narrower but
did not increase the sharp minimum TV on the retained fine/third band.
After projecting the measured child ratios onto the exact original
parent masses, lift the previously feasible dual prices to every child.
The lifted dual has the same objective and remains feasible on every
refined union interval.

Check separately whether the **previous** optimal overlap allocation can
be carried forward without moving it between old union intervals. A
second capacity problem caps all refined descendants by each old overlap
allocation and supplies independent primal and dual witnesses.

* Both forward and inverse number can keep that original allocation,
  with lost overlap below 1.5e−14 in each mesh.
* Forward excess cannot. Under the same old allocation, approximately
  0.00979127 of normalized overlap mass cannot be placed in the children.
  The other refined optimum remains within 3e−14 of the lifted old dual,
  so a different allocation recovers the whole optimum. This does not
  assert uniqueness or that the chosen new solution moves the smallest
  possible mass between old union intervals.

The empirical normalized child counts differ only slightly from the
exact-parent projection; the sum of the two panel TV corrections bounds
any change in minimum compatible TV due to that projection. The result
concerns the minimum compatible TV, not the realized continuous profile.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_cold_native_overlap_inheritance.py --root results \
 --accepted results/early-pole-cold-native-placement-eight-2026-11-28.json
```
