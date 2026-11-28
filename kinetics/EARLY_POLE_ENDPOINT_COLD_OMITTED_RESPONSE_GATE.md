# How much omitted response could change the cold-band conclusion?

The S4/S8/S16 native response profiles are known only on the 17 common display bins inside the retained cold band. Let (P,Q) denote the normalized profiles **conditional on that band**, and (e,f) the fractions of their respective full responses outside it. If the two exterior distributions are unmeasured, their shapes cannot be inferred from the interior. For specified (e,f), the sharp lower and upper possibilities for the full TV distance are
```
L = (sum |(1-e) P - (1-f) Q| + |e-f|)/2
U = (sum |(1-e) P - (1-f) Q| + e+f)/2 .
```
Identical normalized exterior shapes attain the minimum; disjoint exterior shapes attain the maximum. These are extremal mathematical completions, not simulated tails.

The threshold is sharp even if the three unknown fractions differ, provided each is bounded by the **same cap**. For fixed middle retained mass (b=1-f), the lower bound for the preceding comparison can be written as `max(a,b) - sum min(a P_i,b Q_i)`. Increasing (a) up to (b) cannot increase this expression; beyond (b) the overlap can increase no faster than (a), so (a=b) minimizes it. The upper bound for the following comparison is `1 - sum min(b Q_i,c R_i)`, nondecreasing as either omitted fraction grows. Thus the worst case sets all three fractions at the cap, makes the first two exterior shapes identical and the third disjoint. The **smallest cap that permits erasing** the measured decrease is that decrease divided by one plus itself. Above it the full-TV ordering can reverse. These are native profiles, not distances optimized with prescribed moments.

| Response | Cap per outside fraction permitting erasure of S4/S8 → S8/S16 TV contraction | Extra S16 versus S8 exterior rate that could erase the integrated cold decrement, relative to S8 cold rate |
| --- | ---: | ---: |
| Forward number | 0.339669% | 0.0205700% |
| Inverse number | 0.339284% | 0.0241546% |
| Forward excess | 0.327211% | 0.139748% |

The second column is a hypothetical *cap on each fraction*, whereas the last column is a difference of unnormalized rates. The exact interval formula applies separately to arbitrary unequal fractions when those fractions are measured. The thresholds do not show that any such tail exists; they specify why a bound on omitted support is required before extending the cold-band finding to the whole response, let alone to a primordial spectrum.

```sh
cd kinetics
python -W error -m unittest test_early_pole_endpoint_cold_omitted_response_gate
python -W error early_pole_endpoint_cold_omitted_response_gate.py \
  results/early-pole-cold-spatial-population-factor-bounds-2026-11-28.json
```
