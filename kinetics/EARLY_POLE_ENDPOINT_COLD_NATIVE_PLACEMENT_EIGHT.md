# Direct placement refinement on the retained fine and third D2 band

At the same exact archived endpoint states, split at the midpoint the eight
native panels with the largest contribution to the forward-number placement
allowance. Four belong to fine and four to third. On each of two independent
quadrature meshes, evaluate eight new children on each grid. Retain all other
archived rows byte for byte and require the rebuilt whole band to reproduce
the unchanged optical metadata and accepted full rate receipt.

The direct fine/third normalized masses remain free to move within their
new, smaller native subpanels. Compute lower and upper extremal W1 and an
independent minimum compatible TV via common-mass capacities. For the
forward number, inverse number and forward centered excess, compare the
actual result with a second ledger that projects the measured child ratios
onto the exact original masses in their respective parent cells. This
removes parent quadrature drift. The diameter times the sum of the two
normalized TV corrections bounds the possible change in either W1 extremum.
The statically selected "guarantee" is solely the midpoint placement
allowance at fixed original parent masses, not a guarantee about actual rates.

The joined complete pairs have 53 fine and 85 third panels. For forward
number the minimum compatible W1 rises from 6.38274e−8 to 8.77987e−8,
while the sharp maximum falls from 7.02073e−7 to 4.93628e−7. The
forward placement allowance falls by 47.45%, and the ratio of upper
to lower extremal W1 narrows from about 11 to 5.62. Inverse number
and forward excess tighten in the same sense. For all three observables
the minimum compatible TV is unchanged at the reported precision.
The forward mean-q contrast still straddles zero, now approximately
[−4.299e−7,+8.457e−8]. RMS and mean speed at the fixed mass also remain
sign-indeterminate. Replacing empirical child totals with exact parent
conservation bounds the change of either W1 extremum below 9.4e−19 using
the archived finer quadrature, so it does not explain the tightening.
This does not measure physical momentum transfer.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error early_pole_endpoint_cold_native_placement_refinement.py \
 --root results \
 --archive "results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --design "results/early-pole-cold-native-transport-split-targets-$d.json"
```
