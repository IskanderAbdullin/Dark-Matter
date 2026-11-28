# Check numerical sensitivity of the existing conditional fermion floor

At the matched corrected 7.394049 MeV endpoints, repeat the
[piecewise radiation-clock comparison](INDEPENDENT-WARM-SIGNED-FERMION-PIECEWISE-CLOCK-2026-11-26.md).
Use the retained stable fermion count in each exact hashed
snapshot. Do not credit an old exterior or any unknown heavy
number. Transport the original total energy allowance for the
specified lepton and effective charged-pion channels by the
same conservative b³ρ convention, from the 10.442 MeV anchor,
with b=1.4128443 at this later temperature. Retain the previous
subpole and high-loss bounds and their spectral support premises.

| Corrected branch | Fixed-number starting proxy Ωh² | Fine conditional floor at cap 1/2 |
| --- | ---: | ---: |
| q32 / 12 intervals | 75.74274 | 7.929598 |
| q32 / 24 intervals | 75.73549 | 7.929195 |
| q32 / 48 intervals | 75.73176 | 7.928989 |
| q64 / 24 intervals | 75.73547 | 7.929194 |

Varying the time resolution and tested support here changes the
fine comparison by less than 0.008% of its value. The floor is
also evaluated with 350 and 1,000 comparison intervals; both
give weaker values than 2,800. This is a sensitivity of the
**same comparator**, not a statement that any model is excluded
irrespective of its missing reactions. The gross energy ceiling
would need to be roughly 1.877 times the assumed value before
this particular sufficient comparison at cap 1/2 no longer
exceeds .120, with all other premises fixed. That multiplier is
not a limit on unseen sources.

The reconstructed retained light occupation at the anchor is
below 4.11e-5. This snapshot observation does not imply the
required global post-anchor light **and heavy** cap at one-half
per spin. The late full hadronic rates, outside incoming phase,
and independently evolving heavy distribution are not supplied
by this result.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_corrected_anchor_signed_floor_reference.py EOS.dat \
  "$D"/warm-corrected-q32-time12-step12-2026-11-26.json.gz \
  "$D"/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
  "$D"/warm-corrected-q32-time48-step48-2026-11-26.json.gz \
  "$D"/warm-multishell-corrected-q64-stage-24-2026-11-26.json.gz
~~~

[Four matched conditional controls](results/warm-corrected-anchor-signed-floor-2026-11-26.json).
