# High-reheat paired evolution from initially empty dark states

Repeat the 50 MeV-to-12 MeV 48-step physical benchmark using the
same visible source orders, conservative EOS and momentum support,
but insert the signed stable pair channel from the first filled
interval. Retain prompt return, below-pole elasticity, visible
feedback and the separate exterior contribution.

| 48-step visible source mode | Pair conversion | Inclusive Ωχh² | Inclusive Ωφh² |
| --- | --- | ---: | ---: |
| Vacuum visible | No | 88.18921 | 1.45791e-4 |
| Signed visible | No | 88.18859 | 1.45791e-4 |
| Vacuum visible | Yes | 87.18774 | 1.50097e-4 |
| Signed visible | Yes | 87.18726 | 1.50096e-4 |

Integrated stable depletion from this additional channel is only
about 1.14% relative to the matched no-pair run at this endpoint.
The scalar count rises about 2.95%. The gross reverse-to-forward
pair ratio reaches 2.94 during the evolution; a forward-only
depletion estimate would not be the calculation tested here.
The combined energy ledger agrees to better than 6e-14 in
initial visible-energy units. This is a single matched 48-step
mesh at the physical benchmark, not time-step convergence for
all warm channels.

The remaining stable fixed-count proxy is about 726 times .120.
This is not the saturated late relic abundance; the calculation
stops while interactions can still evolve. It does not include
heavy-parent recycling or vacuum scalar 2→4 as collision operators.

```sh
export OPENBLAS_NUM_THREADS=1
common=(python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --visible-backreaction --trace-exterior --with-elastic --empty-pulse 1e-8 --empty-stop .24 --support-low .005 --vacuum-scaled-residual --allow-empty-shells --include-zero-cell --pair-conversion --steps 48)
"${common[@]}" --matched-vacuum-order --save-snapshots warm_vacuum_shape.json > warm_vacuum.json
"${common[@]}" --quantum-visible-source --save-snapshots warm_quantum_shape.json > warm_quantum.json
```
