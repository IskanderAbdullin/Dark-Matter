# Evolve archive-moment-matched packets near the excess-current reversal

The frozen center scan in the prior species-exchange note gives a
near-cut sign bracket but no history. Pick two nearby hard centers
on opposite sides of the .005 GeV half-width bracket, .710 GeV
and .713 GeV. At .713 repeat with the .0025 GeV half-width. All
three use the same original exterior receipt counts and mean
birth energies. The reconstructed daughter energy cutoff
is physically identical, 0.7148213914435242 GeV, and the
finite momentum pad ends at 1.05 GeV.

The archived later phase remains fixed in visible temperature
and is used to initialize a short dark prompt-only test. Advance
the daughter moments with recorded positive implicit
roots, stationary parent closure, no occupation floor and no
hidden subdivision fallback. Verify species number, total
energy, physical W and independently integrated angular current.

For the 5 MeV half-width cases the initial W values are
5.0720451e-18 and 7.1556171e-16 of their initial dark energy.
The initial signed currents at the same physical cut are
+3.2111421e-9 and −6.6398478e-9 per Hubble time, respectively.
These do not imply the signs stay fixed after the source reacts.

| Center / half-width (GeV) | Initial W / rho | Final W / rho | Current signs at Hdt = 0, 3e−8, 1e−7 |
| --- | ---: | ---: | --- |
| .710 / .005 | 5.0720e−18 | 5.6154e−17 | +, +, − |
| .713 / .005 | 7.1556e−16 | 4.5231e−16 | −, −, − |
| .713 / .0025 | 7.7283e−17 | 1.0167e−16 | +, +, − |

The .713 pair has bitwise equal measured initial species
number and energy, with identical H at recorded precision.
The broad packet loses 36.8% of its initial W; the narrow
packet gains 31.6% by the endpoint, despite a negative current
there. The .710 case gains 11.07 times in W relative to
its very small start and its terminal current is only
−3.36e−11 of initial dark energy per Hubble time.
The initial narrow .713 current is +2.54e−9 and its
terminal current is −1.75e−10 in the same units.
Each of the six intervals uses three fixed positive backward
Euler substeps. The recorded change brackets are not
continuum reversal times, nor do they imply that only one
reversal is possible within an interval. In particular the
weak .710 terminal sign should be checked by time refinement.

The change of width is a change in the assumed reconstruction
of the same archive moments, not an estimate of a physical
width. This test supplies neither an unrestricted tail bound
nor evidence for the cosmological production mechanism.
The cutoff is fixed in physical energy and expansion itself
is not advanced in this comparison. See the [frozen
scan](INDEPENDENT-WARM-LATE-TAIL-SPECIES-EXCHANGE-2026-11-26.md)
for the earlier width-dependent initial behavior.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
r() { python -W error warm_late_near_crossing_shape_evolution_reference.py EOS.dat "$a" "$b" "$@"; }
r near710 --hard .710
r near713 --hard .713
r near713narrow --hard .713 --halfwidth .0025
python -W error warm_late_near_crossing_shape_audit_reference.py near710 near713 near713narrow
~~~

[Audit](results/warm-late-near-cut-shapes-audit-2026-11-26.json); compressed complete accepted histories: [710-3](results/warm-late-near-cut-710-3-2026-11-26.json.gz); [713-3](results/warm-late-near-cut-713-3-2026-11-26.json.gz); [713-narrow-3](results/warm-late-near-cut-713-narrow-3-2026-11-26.json.gz).
