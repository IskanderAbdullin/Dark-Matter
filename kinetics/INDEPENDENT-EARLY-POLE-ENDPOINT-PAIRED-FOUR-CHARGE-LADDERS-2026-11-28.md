# Independent fine and third charge ladders share a late downward drift

Repeat the [fine four-charge-subcell history](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-FOUR-CHARGE-SUBCELLS-2026-11-28.md) on the independent third scalar grid. Each ladder holds its scalar grid fixed (67 fine or 97 third bins), starts from empty daughters, and uses the same fixed-time joins. The third four-subcell card has 97 scalar, 1248 daughter and 1408 parent nodes or cells. Both time arms and the paired-number, Pauli and continuity budgets pass for each third continuation segment.

| Charge subcells per archived cell | Fine yield at x=.05 | Third yield at x=.05 | Third relative to fine |
| ---: | ---: | ---: | ---: |
| 1 | 3.78441807e−4 | 3.76522010e−4 | −.5073% |
| 2 | 3.71312973e−4 | 3.65886183e−4 | −1.4615% |
| 4 | 3.60874791e−4 | 3.57755392e−4 | −.8644% |

From two to four subcells both final yields fall, by 2.811% fine and 2.222% third. The two four-subcell charts agree more closely with each other than either agrees with its own coarser calculation. These small peer differences cannot be used to bound the common drift. The third successive decrement ratio is .7645, while the fine ratio is 1.464. This contrast does not support a shared asymptotic power law over the observed levels. The third last fixed-start time drift is 3.74e−5 relative, not a global time error bar.

The third final scalar pole mass .0008621801 is .584% below its original archive, and its normalized scalar histogram TV on the common third grid is .01022. At the upper cold q edge and cosine zero, its moment-only kinematic interval is E/T=21.789–21.802, versus 20.526–20.546 for fine d4. Over all cold q and angles the third maximum threshold is 29.62, below but approaching the archived parent support ceiling near 32. A kinematically available tail is not evidence of appreciable occupancy there. None of these intervals is a physical rate or an optical species-population bound.

Reproduce the third empty start as `early_pole_cell_endpoint_startup_pilot.py third --div 4 --steps 2048 --one` and later restarts with coarse steps 512, 1024, 1024, 1024 at x=.01, .02, .04, .05. The [third chain receipt](results/early-pole-cell-endpoint-third-d4-to-x05-2026-11-28.json), [spectrum audit](results/early-pole-cell-endpoint-third-d4-spectrum-2026-11-28.json), and [threshold receipt](results/early-pole-cell-endpoint-third-d4-threshold-2026-11-28.json) preserve the evidence. Exact restart occupations are archived. To combine the two ladders:
```sh
python -W error early_pole_endpoint_paired_charge_ladders.py \
  results/early-pole-cell-endpoint-{fine,third}-charge-ladder-2026-11-28.json
```
The [paired report](results/early-pole-cell-endpoint-paired-charge-ladders-2026-11-28.json) gives the same comparison at all four fixed times. It is not a combined scalar grid, a continuum extrapolation, or a spectrum ready for CLASS.
