# How much omitted response could change the cold-band conclusion?

The S4/S8/S16 native response profiles are known only on the 17 common display bins inside the retained cold band. Let (P,Q) denote the normalized profiles **conditional on that band**, and (e,f) the fractions of their respective full responses outside it. If the two exterior distributions are unmeasured, their shapes cannot be inferred from the interior. For specified (e,f), the sharp lower and upper possibilities for the full TV distance are
```
L = (sum |(1-e) P - (1-f) Q| + |e-f|)/2
U = (sum |(1-e) P - (1-f) Q| + e+f)/2 .
```
Identical normalized exterior shapes attain the minimum; disjoint exterior shapes attain the maximum. These are extremal mathematical completions, not simulated tails.

If the same exterior fraction is assumed on all three grids, the **smallest possible** fraction permitting the difference between consecutive resolution steps to disappear is the observed cold-only decrease divided by one plus that decrease. To attain equality, give the first two grids identical exterior shapes and the third a disjoint one. Above the threshold the full-TV ordering can reverse. The following values use the measured native profiles, not distances optimized with prescribed moments.

| Response | Common outside fraction that could erase S4/S8 → S8/S16 TV contraction | Extra S16 versus S8 exterior rate that could erase the integrated cold decrement, relative to S8 cold rate |
| --- | ---: | ---: |
| Forward number | 0.339669% | 0.0205700% |
| Inverse number | 0.339284% | 0.0241546% |
| Forward excess | 0.327211% | 0.139748% |

The second column is a hypothetical *common fraction*, whereas the last column is a difference of unnormalized rates. The exact interval formula applies separately to arbitrary unequal fractions when those fractions are measured. The thresholds do not show that any such tail exists; they specify why a bound on omitted support is required before extending the cold-band finding to the whole response, let alone to a primordial spectrum.

```sh
cd kinetics
python -W error -m unittest test_early_pole_endpoint_cold_omitted_response_gate
python -W error early_pole_endpoint_cold_omitted_response_gate.py \
  results/early-pole-cold-spatial-population-factor-bounds-2026-11-28.json
```
